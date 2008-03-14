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

package affymetrix.gcos.exp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

import affymetrix.gcos.FileIO;
import affymetrix.gcos.TagValuePair;

/** Provides parsing and storage for EXP file data. */
public class EXPFileData {

	/** The first line of an EXP file */
	private static final String EXP_HEADER_LINE_1 = "Affymetrix GeneChip Experiment Information";

	/** The second line of an EXP file */
	private static final String EXP_HEADER_LINE_2 = "Version";

	/** The sample information section name. */
	private static final String SAMPLE_SECTION_NAME = "[Sample Info]";

	/** The scanner section name. */
	private static final String SCANNER_SECTION_NAME = "[Scanner]";

	/** The fluidics section name. */
	private static final String FLUIDICS_SECTION_NAME = "[Fluidics]";

	/** The array type tag name. */
	private static final String ARRAY_TYPE_TAG = "Chip Type";

	/** The protocol tag name. */
	private static final String PROTOCOL_TAG = "Protocol";

	/** The station tag name. */
	private static final String STATION_TAG = "Station";

	/** The name of the EXP file */
	private String fileName;

	/** Creates a new instance of EXPFileData */
	public EXPFileData() {
		fileName = "";
		clear();
	}

	/**
	 * Sets the file name.
	 * 
	 * @param name
	 *          The name of the EXP file.
	 */
	public void setFileName(String name) {
		fileName = name;
	}

	/**
	 * Gets the file name.
	 * 
	 * @return The file name.
	 */
	public String getFileName() {
		return fileName;
	}

	/** The array type in the EXP file */
	private String arrayType;

	/**
	 * Gets the array type.
	 * 
	 * @return The probe array type in the mask file.
	 */
	public String getArrayType() {
		return arrayType;
	}

	/** A Vector of scan parameters. */
	private List<TagValuePair> scan = new ArrayList<TagValuePair>();

	/** Gets the number of scan parameters. */
	public int getNumScanParameters() {
		return scan.size();
	}

	/**
	 * Gets the scan parameter for the given index.
	 * 
	 * @param index
	 *          The zero based index to the parameter vector.
	 * @return The tag/value parameter.
	 */
	public TagValuePair getScanParameter(int index) {
		return scan.get(index);
	}

	/** A Vector of hyb parameters. */
	private List<TagValuePair> hyb = new ArrayList<TagValuePair>();

	/** Gets the number of hyb parameters. */
	public int getNumHybParameters() {
		return hyb.size();
	}

	/**
	 * Gets the hyb parameter for the given index.
	 * 
	 * @param index
	 *          The zero based index to the parameter vector.
	 * @return The tag/value parameter.
	 */
	public TagValuePair getHybParameter(int index) {
		return hyb.get(index);
	}

	/** A Vector of sample parameters. */
	private List<TagValuePair> sample = new ArrayList<TagValuePair>();

	/** Gets the number of sample parameters. */
	public int getNumSampleParameters() {
		return sample.size();
	}

	/**
	 * Gets the sample parameter for the given index.
	 * 
	 * @param index
	 *          The zero based index to the parameter vector.
	 * @return The tag/value parameter.
	 */
	public TagValuePair getSampleParameter(int index) {
		return sample.get(index);
	}

	/**
	 * Reads the contents of the file.
	 * 
	 * @return True if successful
	 */
	public boolean read() {
		BufferedReader b = null;
		clear();
		try {
			b = new BufferedReader(new FileReader(fileName));

			// The first two lines are the header.
			String str = FileIO.readNextLine(b);
			if (!str.startsWith(EXP_HEADER_LINE_1)) {
				return false;
			}
			str = FileIO.readNextLine(b);
			if (!str.startsWith(EXP_HEADER_LINE_2)) {
				return false;
			}

			// The possible sections.
			final int NO_SECTION = 0;
			final int SAMPLE_SECTION = 1;
			final int FLUIDICS_SECTION = 2;
			final int SCANNER_SECTION = 3;
			int currentSection = NO_SECTION;

			// The remaining are the sample, fluidics and scanner sections
			boolean captureAll = false;
			while ((str = FileIO.readNextLine(b)) != null) {
				// Check for the start of each section.
				if (str.startsWith(SAMPLE_SECTION_NAME) == true) {
					currentSection = SAMPLE_SECTION;
					continue;
				}

				else if (str.startsWith(FLUIDICS_SECTION_NAME) == true) {
					currentSection = FLUIDICS_SECTION;
					continue;
				}

				else if (str.startsWith(SCANNER_SECTION_NAME) == true) {
					currentSection = SCANNER_SECTION;
					continue;
				}

				// Parse the line into a name/value pair.
				TagValuePair param = new TagValuePair();
				String[] paramStr = str.split("\t");
				param.setTag(paramStr[0]);
				if (paramStr.length == 1) {
					param.setValue("");
				}
				else {
					param.setValue(paramStr[1].trim());
				}

				// Take everything between the "Protocol" and "Station" tags. These
				// may include error messages in the tag name with blank values.
				if ((currentSection == FLUIDICS_SECTION) && (param.getTag().compareTo(PROTOCOL_TAG) == 0)) {
					captureAll = true;
				}
				else if ((currentSection == FLUIDICS_SECTION) && (param.getTag().compareTo(STATION_TAG) == 0)) {
					captureAll = false;
				}

				// If the value is blank then skip it.
				if ((param.getValue().length() == 0) && (captureAll == false)) {
					continue;
				}

				// Double check if the "Protocol" tag is blank and skip it.
				if ((param.getTag().compareTo(PROTOCOL_TAG) == 0) && (param.getValue().length() == 0)) {
					continue;
				}

				// Check for the special array type line
				if ((currentSection == SAMPLE_SECTION) && (param.getTag().compareTo(ARRAY_TYPE_TAG) == 0)) {
					arrayType = param.getValue();
				}

				// Store the parameter into the sample section.
				else if (currentSection == SAMPLE_SECTION) {
					sample.add(param);
				}

				// Store the parameter into the fluidics section.
				else if (currentSection == FLUIDICS_SECTION) {
					hyb.add(param);
				}

				// Store the parameter into the scanner section.
				else if (currentSection == SCANNER_SECTION) {
					scan.add(param);
				}
			}
			return true;
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

	/**
	 * Checks for the existance of a file.
	 * 
	 * @return True if the file exists
	 */
	public boolean exists() {
		return new File(fileName).exists();
	}

	/** Clears memory associated with the class */
	public void clear() {
		arrayType = "";
		scan.clear();
		hyb.clear();
		sample.clear();
	}
}
