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

package affymetrix.fusion.cdf;

import java.io.File;

import affymetrix.gcos.cdf.CDFFileData;
import affymetrix.gcos.cdf.CDFProbeSetInformation;
import affymetrix.gcos.cdf.CDFQCProbeSetInformation;

/** This class provides reading and storage capabilities for the CDF file. */
public class FusionCDFData {

	/** The file header object. */
	private FusionCDFHeader header;

	/** Gets the header. */
	public FusionCDFHeader getHeader() {
		return header;
	}

	/** The GCOS CDF file object. */
	private CDFFileData gcosFile;

	/** The CDF file name (full path). */
	private String fileName;

	/** Gets the file name. */
	public String getFileName() {
		return fileName;
	}

	/** Sets the file name. */
	public void setFileName(final String value) {
		fileName = value;
	}

	/** Gets the error message. */
	public String getError() {
		if (gcosFile != null) {
			return gcosFile.getError();
		}
		return null;
	}

	/**
	 * Gets the name of a probe set.
	 * 
	 * @param index
	 *          The index to the probe set name of interest.
	 * @return The probe set name.
	 */
	public String getProbeSetName(final int index) {
		if (gcosFile != null) {
			return gcosFile.getProbeSetName(index);
		}
		return null;
	}

	/**
	 * Gets the chip type (probe array type) of the CDF file.
	 * 
	 * @return The chip type. This is just the name (without extension) of the CDF file.
	 */
	public String getChipType() {
		if (gcosFile != null) {
			return gcosFile.getChipType();
		}
		return null;
	}

	/** Creates the GCOS objects for reading the files. */
	private void prepareGCOSObjectsForReading() {
		gcosFile = new CDFFileData();
		header = new FusionCDFHeader();
		gcosFile.setFileName(fileName);
	}

	/**
	 * Reads the entire file.
	 * 
	 * @return True if successful.
	 */
	public boolean read() {
		prepareGCOSObjectsForReading();
		final boolean status = gcosFile.read();
		if (status == true) {
			header.setGCOSObject(gcosFile.getHeader());
		}
		return status;
	}

	/**
	 * Reads the header of the file only.
	 * 
	 * @return True if successful.
	 */
	public boolean readHeader() {
		prepareGCOSObjectsForReading();
		final boolean status = gcosFile.readHeader();
		if (status == true) {
			header.setGCOSObject(gcosFile.getHeader());
		}
		return status;
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
		if (gcosFile == null) {
			gcosFile = new CDFFileData();
		}
		gcosFile.setFileName(fileName);
		return gcosFile.isXDACompatibleFile();
	}

	/**
	 * Gets the probe set type for non-qc probe sets.
	 * 
	 * @param index
	 *          The index to the probe set of interest.
	 * @return The type of probe set.
	 */
	public int getProbeSetType(final int index) {
		if (gcosFile != null) {
			return gcosFile.getProbeSetType(index);
		}
		return FusionGeneChipProbeSetType.UnknownProbeSetType;
	}

	/**
	 * Gets the probe set information.
	 * 
	 * @param index
	 *          The index to the probe set of interest.
	 * @param info
	 *          The probe set information.
	 */
	public void getProbeSetInformation(final int index, final FusionCDFProbeSetInformation info) {
		info.clear();
		if (gcosFile != null) {
			final CDFProbeSetInformation probeSet = gcosFile.getProbeSetInformation(index);
			info.setGCOSObject(probeSet);
		}
	}

	/**
	 * Gets the QC probe set information by index.
	 * 
	 * @param index
	 *          The index to the QC probe set of interest.
	 * @param info
	 *          The QC probe set information.
	 */
	public void getQCProbeSetInformation(final int index, final FusionCDFQCProbeSetInformation info) {
		info.clear();
		if (gcosFile != null) {
			final CDFQCProbeSetInformation probeSet = gcosFile.getQCProbeSetInformation(index);
			info.setGCOSObject(probeSet);
		}
	}

	/**
	 * Gets the QC probe set information by type.
	 * 
	 * @param qcType
	 *          The type of QC probe set to retrieve.
	 * @param info
	 *          The QC probe set information.
	 */
	public void getQCProbeSetInformationByType(final int qcType, final FusionCDFQCProbeSetInformation info) {
		info.clear();
		if (gcosFile != null) {
			final CDFQCProbeSetInformation probeSet = gcosFile.getQCProbeSetInformationByType(qcType);
			info.setGCOSObject(probeSet);
		}
	}

	/** Creates a new instance of CDFFileData */
	public FusionCDFData() {
		fileName = "";
		clear();
	}

	/** Clears the members. */
	public void clear() {
		gcosFile = null;
	}

}
