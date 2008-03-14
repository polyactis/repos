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

package affymetrix.fusion.chp;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.data.CHPData;
import affymetrix.calvin.data.CHPExpressionEntry;
import affymetrix.calvin.data.CHPGenotypeEntry;
import affymetrix.calvin.data.CHPUniversalEntry;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parsers.CHPFileReader;
import affymetrix.calvin.parsers.InvalidFileTypeException;
import affymetrix.calvin.parsers.InvalidVersionException;
import affymetrix.calvin.utils.AffymetrixGuidType;
import affymetrix.gcos.chp.CHPFileData;

/** Provides storage for data in a GCOS CHP file. */
public class FusionCHPLegacyData extends FusionCHPData {

	/** The file header object */
	private FusionCHPHeader header;

	/** The GCOS CHP file object. */
	private CHPFileData gcosFile;

	/** The Calvin CHP file object. */
	private CHPData calvinFile;

	/** An error message when failed to read a file. */
	private String errorMessage;

	/** Creates a new instance of FusionCHPData */
	private FusionCHPLegacyData() {
		errorMessage = "";
		clear();
	}

	/**
	 * Accessors to header.
	 * 
	 * @return The header data object
	 */
	public FusionCHPHeader getHeader() {
		return header;
	}

	/**
	 * Returns the expression probe set result
	 * 
	 * @param index
	 *          The index to the result object of interest.
	 * @param entry
	 *          The expression result.
	 */
	public void getExpressionResults(int index, FusionExpressionProbeSetResults entry) throws IOException,
			UnsignedOutOfLimitsException {
		entry.clear();
		if (gcosFile != null) {
			entry.setGCOSObject(gcosFile.getExpressionResults(index));
		}
		else if (calvinFile != null) {
			CHPExpressionEntry calvinEntry = new CHPExpressionEntry();
			calvinFile.getEntry(index, calvinEntry);
			entry.setCalvinObject(calvinEntry);
		}
	}

	/**
	 * Returns the genotyping probe set result
	 * 
	 * @param index
	 *          The index to the result object of interest.
	 * @param entry
	 *          The genotyping result.
	 */
	public void getGenotypingResults(int index, FusionGenotypeProbeSetResults entry) throws IOException,
			UnsignedOutOfLimitsException {
		entry.clear();
		if (gcosFile != null) {
			entry.setGCOSObject(gcosFile.getGenotypingResults(index));
		}
		else if (calvinFile != null) {
			CHPGenotypeEntry calvinEntry = new CHPGenotypeEntry();
			calvinFile.getEntry(index, calvinEntry);
			entry.setCalvinObject(calvinEntry);
		}
	}

	/**
	 * Returns the universal (tag array) probe set result
	 * 
	 * @param index
	 *          The index to the result object of interest.
	 * @param entry
	 *          The universal result.
	 */
	public void getUniversalResults(int index, FusionUniversalProbeSetResults entry) throws UnsignedOutOfLimitsException {
		if (gcosFile != null) {
			entry.setGCOSObject(gcosFile.getUniversalResults(index));
		}
		else if (calvinFile != null) {
			CHPUniversalEntry calvinEntry = new CHPUniversalEntry();
			calvinFile.getEntry(index, calvinEntry);
			entry.setCalvinObject(calvinEntry);
		}
	}

	/**
	 * Returns the resequencing results.
	 * 
	 * @param entry
	 *          The resequencing results.
	 */
	public void getResequencingResults(FusionResequencingResults entry) {
		if (gcosFile != null) {
			entry.setGCOSObject(gcosFile.getResequencingResults());
		}
		else if (calvinFile != null) {
			entry.setCalvinObject(calvinFile);
		}
	}

	/**
	 * Error string when the read functions fail.
	 * 
	 * @return A string message describing a read error
	 */
	public String getError() {
		return errorMessage;
	}

	/**
	 * Reads a calvin file.
	 * 
	 * @return True if successful.
	 */
	private boolean readCalvinFile() {
		CHPFileReader calvinReader = new CHPFileReader();
		try {
			calvinFile = new CHPData();
			calvinReader.setFilename(filename);
			calvinReader.read(calvinFile);
			header = new FusionCHPHeader();
			header.setCalvinObject(calvinFile);
			return true;
		}
		catch (FileNotFoundException e) {
			errorMessage = "File not found";
		}
		catch (InvalidVersionException e) {
			errorMessage = "Invalid file version";
		}
		catch (InvalidFileTypeException e) {
			errorMessage = "Invalid file type";
		}
		catch (IOException e) {
			errorMessage = e.getMessage();
		}
		catch (UnsignedOutOfLimitsException e) {
			errorMessage = e.getDescription();
		}
		clear();
		return false;
	}

	/**
	 * Reads the entire file.
	 * 
	 * @return True if successful.
	 */
	@Override
	protected boolean read() {
		gcosFile = new CHPFileData();
		gcosFile.setFileName(filename);
		if (gcosFile.read()) {
			header = new FusionCHPHeader();
			header.setGCOSHeader(gcosFile.getHeader());
			return true;
		}
		else {
			errorMessage = gcosFile.getError();
		}
		gcosFile = null;
		return readCalvinFile();
	}

	/**
	 * Reads the header of the CHP file
	 * 
	 * @return True if successful
	 */
	@Override
	protected boolean readHeader() {
		gcosFile = new CHPFileData();
		gcosFile.setFileName(filename);
		if (gcosFile.readHeader() == true) {
			header = new FusionCHPHeader();
			header.setGCOSHeader(gcosFile.getHeader());
			return true;
		}
		else {
			errorMessage = gcosFile.getError();
		}
		gcosFile = null;
		return readCalvinFile();
	}

	/**
	 * Get the id of the file (only valid for Command Console "calvin" files)
	 * 
	 * @return The unique file id.
	 */
	@Override
	public AffymetrixGuidType getFileId() {
		if (calvinFile != null) {
			return calvinFile.getGenericData().getFileIdentifier();
		}
		return null;
	}

	/**
	 * Determines if the file specified by the FileName property exists.
	 * 
	 * @return True if the file exists.
	 */
	public boolean exists() {
		return new File(filename).exists();
	}

	/** Deallocates any memory used by the class object */
	public void clear() {
		header = null;
		gcosFile = null;
		calvinFile = null;
	}

	/** A class to register the legacy CHP reader. */
	private static class Reg extends FusionCHPDataReg {
		/** Constructor - register the legacy file type. */
		public Reg() {
			super();
			List<AffymetrixGuidType> ids = new ArrayList<AffymetrixGuidType>();
			AffymetrixGuidType guid;
			guid = new AffymetrixGuidType();
			guid.setGuid(CHPData.CHP_EXPRESSION_ASSAY_TYPE);
			ids.add(guid);
			guid = new AffymetrixGuidType();
			guid.setGuid(CHPData.CHP_RESEQUENCING_ASSAY_TYPE);
			ids.add(guid);
			guid = new AffymetrixGuidType();
			guid.setGuid(CHPData.CHP_GENOTYPING_ASSAY_TYPE);
			ids.add(guid);
			guid = new AffymetrixGuidType();
			guid.setGuid(CHPData.CHP_UNIVERSAL_ASSAY_TYPE);
			ids.add(guid);
			guid = new AffymetrixGuidType();
			ids.add(guid);
			setFileTypeIds(ids);
		}

		/**
		 * Creates a legacy CHP object.
		 * 
		 * @return The legacy CHP object.
		 */
		@Override
		public FusionCHPData makeObject() {
			return new FusionCHPLegacyData();
		}
	};

	/** The one and only registration object. This registers the class as a CHP reader. */
	private static Reg reg = null;

	/** Register the reader with the system. */
	public static void registerReader() {
		if (FusionCHPLegacyData.reg == null) {
			FusionCHPLegacyData.reg = new Reg();
		}
	}

	/**
	 * Converts the type to the legacy CHP type.
	 * 
	 * @param chip
	 *          The pointer to the CHP data object.
	 * @return The legacy CHP data type or NULL if not compatible.
	 */
	public static FusionCHPLegacyData fromBase(FusionCHPData chip) {
		if (chip == null) {
			return null;
		}
		String chipName = chip.getClass().getName();
		String genName = FusionCHPLegacyData.class.getName();
		if (chipName.compareTo(genName) == 0) {
			return (FusionCHPLegacyData)chip;
		}
		return null;
	}
}
