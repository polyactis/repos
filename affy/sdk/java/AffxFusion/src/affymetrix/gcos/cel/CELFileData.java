/////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 2.1 of the License,
// or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
/////////////////////////////////////////////////////////////////

package affymetrix.gcos.cel;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.channels.FileChannel.MapMode;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import affymetrix.gcos.FileIO;
import affymetrix.portability.DataSizes;

/** Provides parsing and data storage for GCOS CEL files. */
public class CELFileData {

	/** Flag to read all of the data in a CEL file. */
	public static final int CEL_ALL = 1;

	/** Flag to read only the data from a CEL file. */
	public static final int CEL_DATA = 2;

	/** Flag to read the outlier and data sections from a CEL file. */
	public static final int CEL_OUTLIER = 4;

	/** Flag to read the mask and data sections from a CEL file. */
	public static final int CEL_MASK = 8;

	/** The magic number for XDA CEL files. */
	public static final int CELL_FILE_MAGIC_NUMBER = 64;

	/** Error string */
	private String strError;

	/**
	 * Gets the error.
	 * 
	 * @return The last error message.
	 */
	public String getError() {
		return strError;
	}

	/** The file name. */
	private String fileName;

	/**
	 * Gets the file name.
	 * 
	 * @return The file name.
	 */
	public String getFileName() {
		return fileName;
	}

	/**
	 * Sets the file name.
	 * 
	 * @param value
	 *          The name of the CEL file to read.
	 */
	public void setFileName(String value) {
		fileName = value;
	}

	/** CEL file header data object. */
	private CELFileHeaderData headerData;

	/**
	 * Gets the header.
	 * 
	 * @return The CEL file header object.
	 */
	public CELFileHeaderData getHeader() {
		return headerData;
	}

	/** The entries for each cell (used for text format). */
	private List<CELFileEntryType> entries = new ArrayList<CELFileEntryType>();

	/** Map for masked cell coordinates. */
	private Map<Integer, Boolean> maskedCells = new HashMap<Integer, Boolean>();

	/** Map for outlier coordinates. */
	private Map<Integer, Boolean> outliers = new HashMap<Integer, Boolean>();

	/** CEL file reading state. */
	// private int nReadState;
	/** Flag to determine if masked cell data should be read */
	private boolean readMaskedCells;

	/** Flag to determine if outlier data should be read */
	private boolean bReadOutliers;

	/** A mapped byte buffer for XDA files. */
	private MappedByteBuffer xdaBuffer;

	/**
	 * Gets the X coordinates from index.
	 * 
	 * @param index
	 *          The 0 based index to the entry array.
	 * @return X coordinate
	 */
	public int indexToX(int index) {
		return index % headerData.getCols();
	}

	/**
	 * Gets the Y coordinates from index.
	 * 
	 * @param index
	 *          The 0 based index to the entry array.
	 * @return Y coordinate
	 */
	public int indexToY(int index) {
		return index / headerData.getCols();
	}

	/**
	 * Maps X/Y coordinates to CEL file index.
	 * 
	 * @param x
	 *          The x coordinate
	 * @param y
	 *          The y coordinate.
	 * @return The index to the entry array.
	 */
	public int xyToIndex(int x, int y) {
		return xyToIndex(x, y, headerData.getRows(), headerData.getCols());
	}

	/**
	 * Maps X/Y coordinates to CEL file index.
	 * 
	 * @param x
	 *          The x coordinate.
	 * @param y
	 *          The y coordinate.
	 * @param r
	 *          The number of rows.
	 * @param c
	 *          The number of columns.
	 * @return The index to the intensity arrays.
	 */
	public static int xyToIndex(int x, int y, int r, int c) {
		return ((y * c) + x);
	}

	/**
	 * Retrieves a CEL file entry.
	 * 
	 * @param index
	 *          The index to the CEL file entries.
	 * @param entry
	 *          The CEL file entry.
	 */
	public void getEntry(int index, CELFileEntryType entry) {
		entry.setIntensity(getIntensity(index));
		entry.setStdv(getStdv(index));
		entry.setPixels(getPixels(index));
	}

	/**
	 * Retrieves a CEL file entry.
	 * 
	 * @param x
	 *          The X coordinate.
	 * @param y
	 *          The Y coordinate.
	 * @param entry
	 *          The CEL file entry.
	 */
	public void getEntry(int x, int y, CELFileEntryType entry) {
		getEntry(xyToIndex(x, y), entry);
	}

	/**
	 * Retrieves a CEL file intensity.
	 * 
	 * @param index
	 *          The index to the CEL file entries.
	 * @return The CEL file intensity.
	 */
	public float getIntensity(int index) {
		assert ((index >= 0) && (index < headerData.getCells()));
		if (xdaBuffer != null) {
			ByteBuffer entry = (ByteBuffer)xdaBuffer.position(index * CELFileEntryType.CEL_FILE_ENTRY_SIZE);
			return entry.getFloat();
		}
		else if (entries.size() > index) {
			CELFileEntryType entry = entries.get(index);
			return entry.getIntensity();
		}
		assert (true);
		return 0.0f;
	}

	/**
	 * Retrieves a CEL file intensity.
	 * 
	 * @param x
	 *          The X coordinate.
	 * @param y
	 *          The Y coordinate.
	 * @return The CEL file intensity.
	 */
	public float getIntensity(int x, int y) {
		return getIntensity(xyToIndex(x, y));
	}

	/**
	 * Retrieves a CEL file stdv value.
	 * 
	 * @param index
	 *          The index to the CEL file entries.
	 * @return The CEL file stdv value.
	 */
	public float getStdv(int index) {
		assert ((index >= 0) && (index < headerData.getCells()));
		if (xdaBuffer != null) {
			ByteBuffer entry = (ByteBuffer)xdaBuffer.position((index * CELFileEntryType.CEL_FILE_ENTRY_SIZE)
					+ DataSizes.FLOAT_SIZE);
			return entry.getFloat();
		}
		else if (entries.size() > index) {
			CELFileEntryType entry = entries.get(index);
			return entry.getStdv();
		}
		assert (true);
		return 0.0f;
	}

	/**
	 * Retrieves a CEL file stdv value.
	 * 
	 * @param x
	 *          The X coordinate.
	 * @param y
	 *          The Y coordinate.
	 * @return The CEL file stdv value.
	 */
	public float getStdv(int x, int y) {
		return getStdv(xyToIndex(x, y));
	}

	/**
	 * Retrieves a CEL file pixel count.
	 * 
	 * @param index
	 *          The index to the CEL file entries.
	 * @return The CEL file pixel count.
	 */
	public short getPixels(int index) {
		assert ((index >= 0) && (index < headerData.getCells()));
		if (xdaBuffer != null) {
			ByteBuffer entry = (ByteBuffer)xdaBuffer.position((index * CELFileEntryType.CEL_FILE_ENTRY_SIZE)
					+ DataSizes.FLOAT_SIZE + DataSizes.FLOAT_SIZE);
			return entry.getShort();
		}
		else if (entries.size() > index) {
			CELFileEntryType entry = entries.get(index);
			return entry.getPixels();
		}
		assert (true);
		return 0;
	}

	/**
	 * Retrieves a CEL file pixel count.
	 * 
	 * @param x
	 *          The X coordinate.
	 * @param y
	 *          The Y coordinate.
	 * @return The CEL file pixel count.
	 */
	public short getPixels(int x, int y) {
		return getPixels(xyToIndex(x, y));
	}

	/**
	 * Retrieves a CEL file mask flag.
	 * 
	 * @param x
	 *          The X coordinate.
	 * @param y
	 *          The Y coordinate.
	 * @return True if the feature is masked.
	 */
	public boolean isMasked(int x, int y) {
		return isMasked(xyToIndex(x, y));
	}

	/**
	 * Retrieves a CEL file mask flag.
	 * 
	 * @param index
	 *          The index to the CEL file entries.
	 * @return True if the feature is masked.
	 */
	public boolean isMasked(int index) {
		Boolean masked = maskedCells.get(new Integer(index));
		if (masked != null) {
			return masked.booleanValue();
		}
		return false;
	}

	/**
	 * Retrieves a CEL file outlier flag.
	 * 
	 * @param x
	 *          The X coordinate.
	 * @param y
	 *          The Y coordinate.
	 * @return True if the feature is an outlier.
	 */
	public boolean isOutlier(int x, int y) {
		return isOutlier(xyToIndex(x, y));
	}

	/**
	 * Retrieves a CEL file outlier flag.
	 * 
	 * @param index
	 *          The index to the CEL file entries.
	 * @return True if the feature is an outlier.
	 */
	public boolean isOutlier(int index) {
		Boolean outlier = outliers.get(new Integer(index));
		if (outlier != null) {
			return outlier.booleanValue();
		}
		return false;
	}

	/**
	 * Checks if the file exists.
	 * 
	 * @return True if the file exists.
	 */
	public boolean exists() {
		return new File(fileName).exists();
	}

	/**
	 * Reads the header of the CEL file.
	 * 
	 * @return True if successful.
	 */
	public boolean readHeader() {
		// Read the header, close if failed.
		if (open(true) == false) {
			clear();
			return false;
		}
		return true;
	}

	/**
	 * Reads the CEL file.
	 * 
	 * @return True if successful.
	 */
	public boolean read() {
		return read(true);
	}

	/**
	 * Reads the CEL file.
	 * 
	 * @param bIncludeMaskAndOutliers
	 *          Flag to indicate if the mask and outlier sections should also be read.
	 * @return True if successful.
	 */
	public boolean read(boolean bIncludeMaskAndOutliers) {
		readMaskedCells = bIncludeMaskAndOutliers;
		bReadOutliers = bIncludeMaskAndOutliers;

		// Open the file
		if (!open(false)) {
			clear();
			return false;
		}
		return true;
	}

	/**
	 * Opens and reads the contents of the CEL file.
	 * 
	 * @param bReadHeaderOnly
	 *          Flag to indicate if the header is only to be read.
	 * @return True if successful.
	 */
	private boolean open(boolean bReadHeaderOnly) {
		if (isXDACompatibleFile()) {
			return readXDACel(bReadHeaderOnly);
		}
		else if (isVersion3CompatibleFile()) {
			return readTextCel(bReadHeaderOnly);
		}
		return false;
	}

	/** Clears the members. */
	public void clear() {
		xdaBuffer = null;
		headerData = null;
		entries.clear();
		outliers.clear();
		maskedCells.clear();
	}

	/**
	 * Checks if the file type is XDA.
	 * 
	 * @return True if XDA type.
	 */
	public boolean isXDACompatibleFile() {
		FileInputStream fis = null;
		try {
			fis = new FileInputStream(fileName);
			int m = FileIO.readInt32(fis);
			return (m == CELFileData.CELL_FILE_MAGIC_NUMBER);
		}
		catch (Throwable t) {
			return false;
		}
		finally {
			try {
				if (fis != null) {
					fis.close();
				}
			}
			catch (Exception e) {
			}
		}
	}

	/**
	 * Checks if the file type is version 3.
	 * 
	 * @return True if version 3 type.
	 */
	public boolean isVersion3CompatibleFile() {
		BufferedReader b = null;
		try {
			b = new BufferedReader(new FileReader(fileName));
			String str = b.readLine();
			b.close();
			if (str.startsWith("[CEL]")) {
				return true;
			}
		}
		catch (Throwable t) {
		}
		finally {
			try {
				if (b != null) {
					b.close();
				}
			}
			catch (Exception e) {
			}
		}
		return false;
	}

	/** Creates a new instance of CELFileData */
	public CELFileData() {
		clear();
	}

	/**
	 * Reads an XDA CEL file.
	 * 
	 * @param bReadHeaderOnly
	 *          Flag to indicate if only the header is to be read.
	 * @return True if successful.
	 */
	private boolean readXDACel(boolean bReadHeaderOnly) {

		FileInputStream fis = null;
		try {
			clear();
			fis = new FileInputStream(fileName);

			// Read the header
			int iHeaderBytes = 0;

			// Read the magic number.
			headerData = new CELFileHeaderData();
			headerData.setMagic(FileIO.readInt32(fis));
			iHeaderBytes += DataSizes.INT_SIZE;

			// Check if new type.
			if (!(headerData.getMagic() == CELFileData.CELL_FILE_MAGIC_NUMBER)) {
				strError = "The file does not appear to be the correct format.";
				return false;
			}

			// Read the version
			headerData.setVersion(FileIO.readInt32(fis));
			iHeaderBytes += DataSizes.INT_SIZE;

			// Read the dimensions of the array
			headerData.setRows(FileIO.readInt32(fis));
			iHeaderBytes += DataSizes.INT_SIZE;
			headerData.setCols(FileIO.readInt32(fis));
			iHeaderBytes += DataSizes.INT_SIZE;
			headerData.setCells(FileIO.readInt32(fis));
			iHeaderBytes += DataSizes.INT_SIZE;

			// Read the other members.
			headerData.setHeader(FileIO.readString(fis));
			iHeaderBytes += DataSizes.INT_SIZE;
			iHeaderBytes += headerData.getHeader().length();
			headerData.setAlg(FileIO.readString(fis));
			iHeaderBytes += DataSizes.INT_SIZE;
			iHeaderBytes += headerData.getAlg().length();
			String str = FileIO.readString(fis);
			headerData.setParameters(str);
			iHeaderBytes += DataSizes.INT_SIZE;
			iHeaderBytes += str.length();
			headerData.setMargin(FileIO.readInt32(fis));
			iHeaderBytes += DataSizes.INT_SIZE;
			headerData.setOutliers(FileIO.readInt32(fis));
			iHeaderBytes += DataSizes.INT_SIZE;
			headerData.setMasked(FileIO.readInt32(fis));
			iHeaderBytes += DataSizes.INT_SIZE;
			/* int ival = */FileIO.readInt32(fis); // ignore this value.
			iHeaderBytes += DataSizes.INT_SIZE;

			// Parse the information from the header.
			headerData.parseChipType();
			headerData.parseDatHeader();
			headerData.parseGrid();

			// Read the remaining data.
			if (bReadHeaderOnly) {
				return true;
			}

			// Get a pointer to the memory map.
			FileChannel fc = fis.getChannel();
			long fileSize = fc.size();
			xdaBuffer = fc.map(MapMode.READ_ONLY, iHeaderBytes, fileSize - iHeaderBytes);
			xdaBuffer.order(ByteOrder.LITTLE_ENDIAN);

			// Read the mask and outlier sections
			int iOffset = headerData.getCells() * (CELFileEntryType.CEL_FILE_ENTRY_SIZE);
			if (readMaskedCells) {
				ByteBuffer entry = (ByteBuffer)xdaBuffer.position(iOffset);
				maskedCells.clear();
				for (int iCell = 0; iCell < headerData.getMasked(); iCell++) {
					maskedCells.put(new Integer(xyToIndex(entry.getShort(), entry.getShort())), new Boolean(true));
				}
			}
			iOffset += (headerData.getMasked() * (DataSizes.SHORT_SIZE + DataSizes.SHORT_SIZE));
			if (bReadOutliers) {
				outliers.clear();
				ByteBuffer entry = (ByteBuffer)xdaBuffer.position(iOffset);
				for (int iCell = 0; iCell < headerData.getOutliers(); iCell++) {
					outliers.put(new Integer(xyToIndex(entry.getShort(), entry.getShort())), new Boolean(true));
				}
			}
			return true;
		}
		catch (Throwable t) {
			strError = t.getMessage();
			return false;
		}
		finally {
			try {
				if (fis != null) {
					fis.close();
				}
			}
			catch (Exception e) {
			}
		}
	}

	/**
	 * Reads an ASCII CEL file.
	 * 
	 * @param bReadHeaderOnly
	 *          Flag to indicate if only the header is to be read.
	 * @return True if successful.
	 */
	private boolean readTextCel(boolean bReadHeaderOnly) {
		BufferedReader b = null;
		try {
			clear();

			// Open the file.
			b = new BufferedReader(new FileReader(fileName));

			// Extract a line of header
			String str = b.readLine();

			// Determine the version number
			headerData = new CELFileHeaderData();
			if (str.startsWith("[CEL]")) {
				headerData.setVersion(3);
			}
			else if (str.equals("COLS/ROWS=")) {
				headerData.setVersion(2);
			}
			else {
				strError = "Unrecognized CEL file format.";
				return false;
			}

			// read and store the header
			if (headerData.getVersion() == 2) {
				headerData.setHeader(str);
				str = b.readLine();
				headerData.setHeader(headerData.getHeader() + "\r\n" + str);
				str = b.readLine();

				// Store rows and columns
				int iCols = 0;
				int iRows = 0;
				// int len = ("COLS/ROWS=").length();
				String[] s = headerData.getHeader().substring(headerData.getHeader().indexOf('=') + 1).split(" ");
				if (s == null) {
					s = headerData.getHeader().substring(headerData.getHeader().indexOf('=') + 1).split("\t");
				}
				iCols = Integer.parseInt(s[0].trim());
				iRows = Integer.parseInt(s[1].trim());
				headerData.setCols(iCols);
				headerData.setRows(iRows);
				headerData.setCells(iCols * iRows);
			}
			else {
				// boolean moreSpace = true;
				// Read past [HEADER] to first line of header data
				while ((str = b.readLine()) != null) {
					if (str.startsWith("[HEADER]")) {
						break;
					}
				}
				int iCols = 0;
				int iRows = 0;
				str = b.readLine();
				String[] s = str.split("=");
				iCols = Integer.parseInt(s[1].trim());
				headerData.setCols(iCols);
				headerData.setHeader(str);
				str = b.readLine();
				s = str.split("=");
				iRows = Integer.parseInt(s[1].trim());
				headerData.setRows(iRows);
				headerData.setHeader(headerData.getHeader() + "\r\n" + str);
				headerData.setCells(iRows * iCols);

				// Now read the rest of the header
				boolean moreHeader = true;
				while (moreHeader) {
					str = b.readLine();

					if (str.startsWith("DatHeader=")) {
						headerData.setDatHeader(str.substring(10));
					}

					if (str.startsWith("Algorithm=")) {
						headerData.setAlg(str.substring(10));
					}

					if (str.startsWith("AlgorithmParameters=")) {
						headerData.setParameters(str.substring(20));
						moreHeader = false;
					}
					headerData.setHeader(headerData.getHeader() + "\r\n" + str);
				}
			}

			// Parse information from the header string.
			headerData.parseChipType();
			headerData.parseMargin();
			headerData.parseGrid();

			// Don't continue if just reading the header.
			if (bReadHeaderOnly) {
				return true;
			}

			// Create memory for Mean data.
			entries.clear();

			// int t_x,t_y,t_pixels;
			// float t_mean,t_stdv;

			// Read v2 CEL files
			if (headerData.getVersion() == 2) {
				// Write the Mean data
				CELFileEntryType entry;
				for (int i = 0; i < headerData.getCells(); i++) {
					String[] s = b.readLine().split("\t");
					entry = new CELFileEntryType();
					entry.setIntensity(Float.parseFloat(s[2].trim()));
					entry.setStdv(Float.parseFloat(s[3].trim()));
					entry.setPixels(Short.parseShort(s[4].trim()));
					entries.add(entry);
				}
			}
			else {
				// Advance to the beginning of the Mean data
				while ((str = b.readLine()) != null) {
					if (str.startsWith("[INTENSITY]")) {
						break;
					}
				}
				str = b.readLine();
				str = b.readLine();

				// Read the Mean data
				CELFileEntryType entry;
				for (int i = 0; i < headerData.getCells(); i++) {
					String[] s = b.readLine().split("\t");
					entry = new CELFileEntryType();
					entry.setIntensity(Float.parseFloat(s[2].trim()));
					entry.setStdv(Float.parseFloat(s[3].trim()));
					entry.setPixels(Short.parseShort(s[4].trim()));
					entries.add(entry);
				}

				// Advance to the Masked data
				while ((str = b.readLine()) != null) {
					if (str.startsWith("[MASKS]")) {
						break;
					}
				}
				if (str == null) {
					return true;
				}

				// Read the masked data
				if (readMaskedCells) {
					// Read number of masked cells
					str = b.readLine();
					int nMasked = Integer.parseInt(str.substring(str.indexOf('=') + 1).trim());
					maskedCells.clear();
					headerData.setMasked(nMasked);
					str = b.readLine(); // skip over the header
					for (int i = 0; i < nMasked; i++) {
						String[] s = b.readLine().split("\t");
						maskedCells.put(new Integer(xyToIndex(Integer.parseInt(s[0].trim()), Integer.parseInt(s[1].trim()))),
								new Boolean(true));
					}
				}
				else {
					headerData.setMasked(0);
				}
				// Advance to the Masked data
				while ((str = b.readLine()) != null) {
					if (str.startsWith("[OUTLIERS]")) {
						break;
					}
				}
				if (str == null) {
					return true;
				}

				// Read the outlier data
				if (bReadOutliers) {
					// Read number of outlier cells
					str = b.readLine();
					int nOutliers = Integer.parseInt(str.substring(str.indexOf('=') + 1).trim());
					outliers.clear();
					headerData.setOutliers(nOutliers);
					str = b.readLine(); // skip over the header
					for (int i = 0; i < nOutliers; i++) {
						String[] s = b.readLine().split("\t");
						outliers.put(new Integer(xyToIndex(Integer.parseInt(s[0].trim()), Integer.parseInt(s[1].trim()))),
								new Boolean(true));
					}
				}
				else {
					headerData.setOutliers(0);
				}
			}
			return true;
		}
		catch (Throwable t) {
			clear();
			strError = t.getMessage();
			return false;
		}
		finally {
			try {
				if (b != null) {
					b.close();
				}
			}
			catch (Exception e) {
			}
		}
	}
}
