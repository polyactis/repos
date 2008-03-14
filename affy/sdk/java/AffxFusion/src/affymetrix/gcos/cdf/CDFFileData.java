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

package affymetrix.gcos.cdf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.channels.FileChannel.MapMode;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import affymetrix.gcos.FileIO;
import affymetrix.portability.DataSizes;

/** This class provides reading and storage capabilities for the CDF file. */
public class CDFFileData {

	/** The magic number for an XDA CDF file. */
	private static final int CDF_FILE_MAGIC_NUMBER = 67;

	/** The version number for an XDA CDF file. */
	private static final int CDF_FILE_VERSION_NUMBER = 1;

	private static final String EQ = "=";

	/** The file header object. */
	private CDFFileHeader header;

	/** The list of probe set names. */
	private CDFProbeSetNames probeSetNames;

	/** An array of probe sets. */
	private Vector<CDFProbeSetInformation> probeSets = new Vector<CDFProbeSetInformation>();

	/** An array of QC probe sets. */
	private Vector<CDFQCProbeSetInformation> qcProbeSets = new Vector<CDFQCProbeSetInformation>();

	/** The CDF file name (full path). */
	private String fileName;

	/** Gets the header. */
	public CDFFileHeader getHeader() {
		return header;
	}

	/** Gets the file name. */
	public String getFileName() {
		return fileName;
	}

	/** Sets the file name. */
	public void setFileName(String value) {
		fileName = value;
	}

	/** A string to hold an error message upon read failures. */
	private String strError;

	/** Gets the error message. */
	public String getError() {
		return strError;
	}

	/** Array of file positions for the probe set data. */
	private List<Integer> probeSetPositions = new ArrayList<Integer>();

	/** Array of file positions for the QC probe set data. */
	private List<Integer> qcProbeSetPositions = new ArrayList<Integer>();

	/** A mapped byte buffer for XDA files. */
	private MappedByteBuffer xdaBuffer;

	/**
	 * Opens the file for reading.
	 * 
	 * @param bReadHeaderOnly
	 *          Flag indicating if only the header should be read.
	 * @return True if successful.
	 */
	private boolean open(boolean bReadHeaderOnly) {
		clear();
		if (isXDACompatibleFile()) {
			return readXDAFormat(bReadHeaderOnly);
		}
		else {
			return readTextFormat(bReadHeaderOnly);
		}
	}

	/**
	 * Reads the header from a text CDF file.
	 * 
	 * @param b
	 *          The file buffer reader.
	 * @return True if successful.
	 */
	private boolean readTextHeader(BufferedReader b) {
		String str;
		String[] sarray;
		final String CDFVERSION1 = "GC1.0";
		final String CDFVERSION2 = "GC2.0";
		final String CDFVERSION3 = "GC3.0";

		// Get the CDF section.
		str = FileIO.readNextLine(b);
		if (str.startsWith("[CDF]") == false) {
			strError = "Unknown file format.";
			return false;
		}

		// Get the version number.
		header = new CDFFileHeader();
		sarray = FileIO.readNextLine(b).split(EQ);
		if (sarray[1].startsWith(CDFVERSION1) == true) {
			header.setVersion(1);
		}
		else if (sarray[1].startsWith(CDFVERSION2) == true) {
			header.setVersion(2);
		}
		else if (sarray[1].startsWith(CDFVERSION3) == true) {
			header.setVersion(3);
		}

		// Get the next section.
		str = FileIO.readNextLine(b); // [Chip]
		str = FileIO.readNextLine(b); // name
		sarray = FileIO.readNextLine(b).split(EQ); // rows
		header.setRows(Integer.parseInt(sarray[1].trim()));
		sarray = FileIO.readNextLine(b).split(EQ); // cols
		header.setCols(Integer.parseInt(sarray[1].trim()));
		sarray = FileIO.readNextLine(b).split(EQ); // #ProbeSets
		header.setNumProbeSets(Integer.parseInt(sarray[1].trim()));
		str = FileIO.readNextLine(b); // max ProbeSet number
		header.setNumQCProbeSets(0);
		if (header.getVersion() > 1) {
			sarray = FileIO.readNextLine(b).split(EQ); // #qc ProbeSets
			header.setNumQCProbeSets(Integer.parseInt(sarray[1].trim()));
			sarray = FileIO.readNextLine(b).split(EQ); // The reference string.
			if (sarray.length == 2) {
				header.setReference(sarray[1]);
			}
		}
		return true;
	}

	/**
	 * Read the QC section from a text CDF file.
	 * 
	 * @param b
	 *          The file buffer reader.
	 */
	private void readTextQC(BufferedReader b) {
		// String str;
		String[] sarray;

		// Allocate for the QCProbeSets.
		qcProbeSets.clear();

		// Read the QC probe sets
		for (int iQCProbeSet = 0; iQCProbeSet < header.getNumQCProbeSets(); iQCProbeSet++) {
			CDFQCProbeSetInformation qcProbeSet = new CDFQCProbeSetInformation();

			/* str = */FileIO.readNextLine(b); // label [QCUnit...]
			sarray = FileIO.readNextLine(b).split(EQ); // type
			qcProbeSet.setQCProbeSetType(Short.parseShort(sarray[1].trim()));
			sarray = FileIO.readNextLine(b).split(EQ); // #cells
			qcProbeSet.setNumCells(Integer.parseInt(sarray[1].trim()));
			/* str = */FileIO.readNextLine(b); // cell header

			// Read the QC cells.
			int xqc;
			int yqc;
			byte plenqc;
			List<CDFQCProbeInformation> qccells = new ArrayList<CDFQCProbeInformation>();
			for (int iqccell = 0; iqccell < qcProbeSet.getNumCells(); iqccell++) {
				CDFQCProbeInformation qcCell = new CDFQCProbeInformation();

				sarray = FileIO.readNextLine(b).split(EQ);
				sarray = sarray[1].split("\t");
				xqc = Integer.parseInt(sarray[0].trim());
				yqc = Integer.parseInt(sarray[1].trim());
				plenqc = Byte.parseByte(sarray[3].trim());

				qcCell.setX(xqc);
				qcCell.setY(yqc);
				qcCell.setProbeLength(plenqc);
				qcCell.setBackground(false);
				qcCell.setPMProbe(false);

				qccells.add(qcCell);
			}
			qcProbeSet.setCells(qccells);
			qcProbeSets.add(qcProbeSet);
		}
	}

	/**
	 * Reads a text format CDF file.
	 * 
	 * @param bReadHeaderOnly
	 *          Flag indicating if only the header should be read.
	 * @return True if successful.
	 */
	private boolean readTextFormat(boolean bReadHeaderOnly) {
		BufferedReader b = null;
		try {
			b = new BufferedReader(new FileReader(fileName));

			// Read the header
			if (!readTextHeader(b)) {
				return false;
			}

			// Stop if just reading the header.
			if (bReadHeaderOnly) {
				return true;
			}
			// Read the remaing sections.
			readTextQC(b);
			readTextProbeSets(b);
		}
		catch (Throwable t) {
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
		return true;
	}

	/**
	 * Read the probe sets from the text CDF file.
	 * 
	 * @param b
	 *          The file buffer reader.
	 */
	private void readTextProbeSets(BufferedReader b) {
		String str;
		String[] sarray;

		// Allocate for the probe set names.
		probeSetNames = new CDFProbeSetNames();
		probeSetNames.setSize(header.getNumProbeSets());

		// Allocate for the ProbeSets.
		int iProbeSet = 0;
		probeSets.clear();
		probeSets.setSize(header.getNumProbeSets());

		// Skip until the ProbeSet section is found
		boolean nextProbeSet = true;
		while (nextProbeSet == true) {
			while (true) {
				str = FileIO.readNextLine(b);
				if (str == null) {
					return;
				}
				if ((str.length() > 5) && (str.startsWith("[Unit") == true) && (str.indexOf("Block") == -1)) {
					break;
				}
			}

			boolean expectMisMatch = false;
			CDFProbeSetInformation probeSet = new CDFProbeSetInformation();

			// ProbeSet info.
			probeSet.setIndex(iProbeSet);
			sarray = FileIO.readNextLine(b).split(EQ); // name
			probeSetNames.setName(iProbeSet, sarray[1]);
			sarray = FileIO.readNextLine(b).split(EQ); // direction
			probeSet.setDirection(Byte.parseByte(sarray[1].trim()));
			sarray = FileIO.readNextLine(b).split(EQ); // # Lists
			sarray = sarray[1].split("\t");
			probeSet.setNumCellsPerList((byte)0);
			probeSet.setNumLists(Integer.parseInt(sarray[0].trim()));
			if (sarray.length == 2) {
				probeSet.setNumCellsPerList(Byte.parseByte(sarray[1].trim()));
			}
			sarray = FileIO.readNextLine(b).split(EQ); // # cells
			probeSet.setNumCells(Integer.parseInt(sarray[1].trim()));
			sarray = FileIO.readNextLine(b).split(EQ); // ProbeSet number
			probeSet.setProbeSetNumber(Integer.parseInt(sarray[1].trim()));
			sarray = FileIO.readNextLine(b).split(EQ); // type
			int ival = Integer.parseInt(sarray[1].trim());

			// final int UNKNOWN_TILE = 0;
			final int STANDARD_TILE = 1;
			final int BLOCK_TILE = 2;
			final int GENE_EXPRESSION_TILE = 3;
			// final int CONTROL_TILE = 4;
			final int STANDARD_ALTERNATE_TILE = 5;
			final int STANDARD_VARIANT_TILE = 6;
			final int UNIVERSAL_TILE = 7;
			final int COPY_NUMBER_TILE = 8;
			final int GENOTYPE_CONTROL_TILE = 9;
			final int EXPRESSION_CONTROL_TILE = 10;

			switch (ival) {
			case STANDARD_TILE:
			case STANDARD_ALTERNATE_TILE:
			case STANDARD_VARIANT_TILE:
				probeSet.setProbeSetType(GeneChipProbeSetType.ResequencingProbeSetType);
				break;

			case BLOCK_TILE:
				probeSet.setProbeSetType(GeneChipProbeSetType.GenotypingProbeSetType);
				break;

			case GENE_EXPRESSION_TILE:
				probeSet.setProbeSetType(GeneChipProbeSetType.ExpressionProbeSetType);
				break;

			case UNIVERSAL_TILE:
				probeSet.setProbeSetType(GeneChipProbeSetType.TagProbeSetType);
				break;

			case COPY_NUMBER_TILE:
				probeSet.setProbeSetType(GeneChipProbeSetType.CopyNumberProbeSetType);
				break;

			case GENOTYPE_CONTROL_TILE:
				probeSet.setProbeSetType(GeneChipProbeSetType.GenotypeControlProbeSetType);
				break;

			case EXPRESSION_CONTROL_TILE:
				probeSet.setProbeSetType(GeneChipProbeSetType.ExpressionControlProbeSetType);
				break;

			default:
				probeSet.setProbeSetType(GeneChipProbeSetType.UnknownProbeSetType);
				break;
			}

			sarray = FileIO.readNextLine(b).split(EQ); // # blocks
			probeSet.setNumGroups(Integer.parseInt(sarray[1].trim()));

			// Determine the number of cells per List if not specified
			// in the CDF file.
			if (probeSet.getNumCellsPerList() == 0) {
				if ((probeSet.getProbeSetType() == GeneChipProbeSetType.GenotypingProbeSetType)
						|| (probeSet.getProbeSetType() == GeneChipProbeSetType.ResequencingProbeSetType)
						|| (probeSet.getProbeSetType() == GeneChipProbeSetType.TagProbeSetType)
						|| (probeSet.getProbeSetType() == GeneChipProbeSetType.UnknownProbeSetType)) {
					probeSet.setNumCellsPerList((byte)4);
				}
				else if ((probeSet.getProbeSetType() == GeneChipProbeSetType.ExpressionProbeSetType)
						|| (probeSet.getProbeSetType() == GeneChipProbeSetType.CopyNumberProbeSetType)
						|| (probeSet.getProbeSetType() == GeneChipProbeSetType.GenotypeControlProbeSetType)
						|| (probeSet.getProbeSetType() == GeneChipProbeSetType.ExpressionControlProbeSetType)) {
					if ((probeSet.getNumLists() != 0) && (probeSet.getNumCells() / probeSet.getNumLists() < 255)) {
						probeSet.setNumCellsPerList((byte)(probeSet.getNumCells() / probeSet.getNumLists()));
					}
					else {
						probeSet.setNumCellsPerList((byte)1);
					}
				}
				else {
					probeSet.setNumCellsPerList((byte)1);
				}
			}

			// If this is an expression probe set and we have 2 cells per list set expectMisMatch flag.
			if ((probeSet.getProbeSetType() == GeneChipProbeSetType.ExpressionProbeSetType)
					&& (probeSet.getNumCellsPerList() == 2)) {
				expectMisMatch = true;
			}

			// Get the mutation type if block tile. ignore.
			if ((probeSet.getProbeSetType() == GeneChipProbeSetType.GenotypingProbeSetType) && (header.getVersion() > 1)) {
				str = FileIO.readNextLine(b);
			}

			// Read the blocks.
			Vector<CDFProbeGroupInformation> groups = new Vector<CDFProbeGroupInformation>(probeSet.getNumGroups());
			groups.setSize(probeSet.getNumGroups());
			for (int iGroup = 0; iGroup < probeSet.getNumGroups(); iGroup++) {
				CDFProbeGroupInformation blk = new CDFProbeGroupInformation();
				blk.setGroupIndex(iGroup);
				blk.setProbeSetIndex(iProbeSet);

				str = FileIO.readNextLine(b); // section name - ignore
				sarray = FileIO.readNextLine(b).split(EQ); // name
				blk.setName(sarray[1]);

				if (probeSet.getProbeSetType() == GeneChipProbeSetType.ExpressionProbeSetType) {
					probeSetNames.setName(iProbeSet, sarray[1]);
				}

				str = FileIO.readNextLine(b); // block number - ignore.
				sarray = FileIO.readNextLine(b).split(EQ); // number of Lists.
				blk.setNumLists(Integer.parseInt(sarray[1].trim()));
				sarray = FileIO.readNextLine(b).split(EQ); // number of cells
				blk.setNumCells(Integer.parseInt(sarray[1].trim()));
				sarray = FileIO.readNextLine(b).split(EQ); // start position.
				blk.setStart(Integer.parseInt(sarray[1].trim()));
				sarray = FileIO.readNextLine(b).split(EQ); // stop position
				blk.setStop(Integer.parseInt(sarray[1].trim()));
				blk.setNumCellsPerList(probeSet.getNumCellsPerList());
				if ((probeSet.getProbeSetType() == GeneChipProbeSetType.GenotypingProbeSetType) && (header.getVersion() > 2)) {
					sarray = FileIO.readNextLine(b).split(EQ);
					blk.setDirection(Byte.parseByte(sarray[1].trim()));
				}
				else {
					blk.setDirection(probeSet.getDirection());
				}

				// Read the cells.
				str = FileIO.readNextLine(b); // header
				Vector<CDFProbeInformation> cells = new Vector<CDFProbeInformation>(blk.getNumCells());
				cells.setSize(blk.getNumCells());
				int cellIndex = 0;
				for (int iCell = 0; iCell < blk.getNumCells(); iCell++) {
					CDFProbeInformation cell = new CDFProbeInformation();
					sarray = FileIO.readNextLine(b).split(EQ);
					sarray = sarray[1].split("\t");
					cell.setX(Integer.parseInt(sarray[0].trim()));
					cell.setY(Integer.parseInt(sarray[1].trim()));
					cell.setExpos(Integer.parseInt(sarray[5].trim()));
					cell.setPBase(sarray[8].charAt(0));
					cell.setTBase(sarray[9].charAt(0));
					cell.setListIndex(Integer.parseInt(sarray[10].trim()));

					if (probeSet.getProbeSetType() == GeneChipProbeSetType.ExpressionProbeSetType) {
						cellIndex = (iCell / probeSet.getNumCellsPerList()) * probeSet.getNumCellsPerList();
						if (expectMisMatch && (cell.getPBase() == cell.getTBase())) {
							++cellIndex;
						}
					}
					else {
						cellIndex = (iCell / probeSet.getNumCellsPerList()) * probeSet.getNumCellsPerList();
						cellIndex += (probeSet.getNumCellsPerList() - (iCell % probeSet.getNumCellsPerList()) - 1);
					}

					cells.set(cellIndex, cell);

					if (iCell == 0) {
						blk.setStart(cell.getListIndex());
					}
					else if (iCell == blk.getNumCells() - 1) {
						blk.setStop(cell.getListIndex());
					}
				}
				blk.setCells(cells);
				groups.set(iGroup, blk);
			}
			probeSet.setGroups(groups);
			probeSets.set(iProbeSet, probeSet);
			++iProbeSet;
		}
		return;
	}

	/**
	 * Reads an XDA format CDF file.
	 * 
	 * @param bReadHeaderOnly
	 *          Flag indicating if only the header should be read.
	 * @return True if successful.
	 */
	private boolean readXDAFormat(boolean bReadHeaderOnly) {

		// Open the file.
		FileInputStream fis = null;
		try {
			fis = new FileInputStream(fileName);
		}
		catch (Throwable t) {
			strError = t.getMessage();
			return false;
		}

		// Read the header.
		if (readXDAHeader(fis) == false) {
			return false;
		}

		// Stop if just reading the header.
		if (bReadHeaderOnly) {
			return true;
		}

		// Get a pointer to the memory map.
		FileChannel fc = fis.getChannel();
		int headerOffset = 0;
		try {
			headerOffset = (int)fc.position();
			long fileSize = fc.size();
			xdaBuffer = fc.map(MapMode.READ_ONLY, headerOffset, fileSize - headerOffset);
			xdaBuffer.order(ByteOrder.LITTLE_ENDIAN);
			fis.close();
		}
		catch (Throwable t) {
			strError = t.getMessage();
			return false;
		}

		// Now that the file is mapped, set the file pointers of the members
		int dataOffset = 0;
		probeSetNames = new CDFProbeSetNames();
		probeSetNames.setMap(xdaBuffer, dataOffset);

		// Skip over the probe set names
		dataOffset += (CDFProbeSetNames.MAX_PROBE_SET_NAME_LENGTH * header.getNumProbeSets());

		// Read the qc probe set indicies
		qcProbeSetPositions.clear();
		for (int iqcset = 0; iqcset < header.getNumQCProbeSets(); iqcset++) {
			qcProbeSetPositions.add(iqcset, FileIO.getInt32(xdaBuffer, dataOffset) - headerOffset);
			dataOffset += DataSizes.INT_SIZE;
		}

		// Read the probe set indicies.
		probeSetPositions.clear();
		for (int iset = 0; iset < header.getNumProbeSets(); iset++) {
			probeSetPositions.add(iset, FileIO.getInt32(xdaBuffer, dataOffset) - headerOffset);
			dataOffset += DataSizes.INT_SIZE;
		}
		return true;
	}

	/**
	 * Reads the header of an XDA format CDF file.
	 * 
	 * @param instr
	 *          The file stream object.
	 * @return True if successful.
	 */
	private boolean readXDAHeader(FileInputStream fis) {
		// Extact the magic and version numbers.
		header = new CDFFileHeader();
		header.setMagic(FileIO.readInt32(fis));
		header.setVersion(FileIO.readInt32(fis));

		// Check the values for the right format file.
		if ((header.getMagic() != CDF_FILE_MAGIC_NUMBER) || (header.getVersion() > CDF_FILE_VERSION_NUMBER)) {
			strError = "The file does not appear to be the correct format.";
			return false;
		}

		// Read the remaining header.
		header.setCols(FileIO.readInt16(fis));
		header.setRows(FileIO.readInt16(fis));
		header.setNumProbeSets(FileIO.readInt32(fis));
		header.setNumQCProbeSets(FileIO.readInt32(fis));
		header.setReference(FileIO.readString(fis));

		return true;
	}

	/**
	 * Gets the name of a probe set.
	 * 
	 * @param index
	 *          The index to the probe set name of interest.
	 * @return The probe set name.
	 */
	public String getProbeSetName(int index) {
		return probeSetNames.getName(index);
	}

	/**
	 * Gets the chip type (probe array type) of the CDF file.
	 * 
	 * @return The chip type. This is just the name (without extension) of the CDF file.
	 */
	public String getChipType() {
		if (fileName.length() > 0) {
			int start = fileName.lastIndexOf('\\');
			if (start == -1) {
				start = fileName.lastIndexOf('/');
			}
			int end = fileName.lastIndexOf('.');
			return fileName.substring(start + 1, end);
		}
		return "";
	}

	/**
	 * Reads the entire file.
	 * 
	 * @return True if successful.
	 */
	public boolean read() {
		if (open(false) == false) {
			clear();
			return false;
		}
		return true;
	}

	/**
	 * Reads the header of the file only.
	 * 
	 * @return True if successful.
	 */
	public boolean readHeader() {
		if (open(true) == false) {
			clear();
			return false;
		}
		return true;
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
	 * Determines if a CDF file is of the XDA (binary) format.
	 * 
	 * @return True if XDA format.
	 */
	public boolean isXDACompatibleFile() {
		FileInputStream fis = null;
		try {
			fis = new FileInputStream(fileName);
			int m = FileIO.readInt32(fis);
			return (m == CDFFileData.CDF_FILE_MAGIC_NUMBER);
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
	 * Gets the probe set type for non-qc probe sets.
	 * 
	 * @param index
	 *          The index to the probe set of interest.
	 * @return The type of probe set.
	 */
	public int getProbeSetType(int index) {
		if (probeSets.size() > 0) {
			return probeSets.get(index).getProbeSetType();
		}
		else if (xdaBuffer != null) {
			// The type is the first item in the probe set object.
			int offset = probeSetPositions.get(index);
			return FileIO.getUInt16(xdaBuffer, offset);
		}
		return GeneChipProbeSetType.UnknownProbeSetType;
	}

	/**
	 * Gets the probe set information.
	 * 
	 * @param index
	 *          The index to the probe set of interest.
	 * @return The probe set information.
	 */
	public CDFProbeSetInformation getProbeSetInformation(int index) {
		if (probeSets.size() > 0) {
			return probeSets.get(index);
		}
		else if (xdaBuffer != null) {
			CDFProbeSetInformation probeSet = new CDFProbeSetInformation();
			long dataOffset = probeSetPositions.get(index);
			probeSet.setMap(xdaBuffer, dataOffset, index);
			return probeSet;
		}
		return null;
	}

	/**
	 * Gets the QC probe set information by index.
	 * 
	 * @param index
	 *          The index to the QC probe set of interest.
	 * @return The QC probe set information.
	 */
	public CDFQCProbeSetInformation getQCProbeSetInformation(int index) {
		if (qcProbeSets.size() > 0) {
			return qcProbeSets.get(index);
		}
		else if (xdaBuffer != null) {
			CDFQCProbeSetInformation qcProbeSet = new CDFQCProbeSetInformation();
			long dataOffset = qcProbeSetPositions.get(index);
			qcProbeSet.setMap(xdaBuffer, dataOffset);
			return qcProbeSet;
		}
		return null;
	}

	/**
	 * Gets the QC probe set information by type.
	 * 
	 * @param qcType
	 *          The type of QC probe set to retrieve.
	 * @return The QC probe set information.
	 */
	public CDFQCProbeSetInformation getQCProbeSetInformationByType(int qcType) {
		for (int i = 0; i < header.getNumQCProbeSets(); i++) {
			CDFQCProbeSetInformation info = getQCProbeSetInformation(i);
			if (info.getQCProbeSetType() == qcType) {
				return info;
			}
		}
		return null;
	}

	/** Creates a new instance of CDFFileData */
	public CDFFileData() {
		fileName = "";
		strError = "";
		clear();
	}

	/** Clears the members. */
	public void clear() {
		xdaBuffer = null;
		probeSetPositions.clear();
		qcProbeSetPositions.clear();
		header = null;
		probeSetNames = null;
		probeSets.clear();
		qcProbeSets.clear();
	}

}
