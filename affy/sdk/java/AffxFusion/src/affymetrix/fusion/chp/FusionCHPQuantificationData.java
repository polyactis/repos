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

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.data.CHPQuantificationData;
import affymetrix.calvin.data.GenericData;
import affymetrix.calvin.data.ProbeSetQuantificationData;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.parsers.CHPQuantificationFileReader;
import affymetrix.calvin.utils.AffymetrixGuidType;

/** This class provides storage and reading capabilities for Quantification CHP files */
public class FusionCHPQuantificationData extends FusionCHPData {
	/** Constructor */
	private FusionCHPQuantificationData() {
		chpData = null;
	}

	/** The CHP object. */
	private CHPQuantificationData chpData;

	/**
	 * Get the id of the file (only valid for Command Console "calvin" files)
	 * 
	 * @return The unique file id.
	 */
	@Override
	public AffymetrixGuidType getFileId() {
		return chpData.getGenericData().getFileIdentifier();
	}

	/** Returns the GenericData object associated with a Calvin file, NULL for GCOS files. */
	public GenericData getGenericData() {
		return chpData.getGenericData();
	}

	/**
	 * Gets the name of the algorithm.
	 * 
	 * @return The algorithm name.
	 */
	public String getAlgName() {
		return chpData.getAlgName();
	}

	/**
	 * Gets the algorithm version.
	 * 
	 * @return The version.
	 */
	public String getAlgVersion() {
		return chpData.getAlgVersion();
	}

	/** Sets the array type */
	public String getArrayType() {
		return chpData.getArrayType();
	}

	/**
	 * Gets the algorithm parameters
	 * 
	 * @return The algoirhtm parameters.
	 */
	public List<ParameterNameValue> getAlgParams() {
		return chpData.getAlgParams();
	}

	/**
	 * Gets the summary parameters
	 * 
	 * @return The summary parameters.
	 */
	public List<ParameterNameValue> getSummaryParams() {
		return chpData.getSummaryParams();
	}

	/** Gets the number of entries (probe sets) */
	public int getEntryCount() {
		return chpData.getEntryCount();
	}

	/**
	 * Gets the Quantification/probe set name data for the given row.
	 * 
	 * @param index
	 *          The row index.
	 * @return The entry.
	 */
	public ProbeSetQuantificationData getQuantificationEntry(int index) throws IOException, UnsignedOutOfLimitsException {
		return chpData.getQuantificationEntry(index);
	}

	/**
	 * Reads the CHP file.
	 * 
	 * @return True if successful.
	 */
	@Override
	protected boolean read() {
		CHPQuantificationFileReader reader = new CHPQuantificationFileReader();
		chpData = new CHPQuantificationData();
		reader.setFilename(filename);
		try {
			reader.read(chpData);
			return true;
		}
		catch (Throwable t) {
			return false;
		}
	}

	/**
	 * Reads the header of the CHP file
	 * 
	 * @return True if successful
	 */
	@Override
	protected boolean readHeader() {
		return read();
	}

	/** A class to register the Quantification CHP reader. */
	private static class Reg extends FusionCHPDataReg {
		/** Constructor - register the Quantification file type. */
		public Reg() {
			super();
			List<AffymetrixGuidType> ids = new ArrayList<AffymetrixGuidType>();
			AffymetrixGuidType guid = new AffymetrixGuidType();
			guid.setGuid(CHPQuantificationData.CHP_QUANTIFICATION_TYPE);
			ids.add(guid);
			setFileTypeIds(ids);
		}

		/**
		 * Creates a Quantification CHP object.
		 * 
		 * @return The Quantification CHP object.
		 */
		@Override
		public FusionCHPData makeObject() {
			return new FusionCHPQuantificationData();
		}
	};

	/** The one and only registration object. This registers the class as a CHP reader. */
	private static Reg reg = null;

	/** Register the reader with the system. */
	public static void registerReader() {
		if (FusionCHPQuantificationData.reg == null) {
			FusionCHPQuantificationData.reg = new Reg();
		}
	}

	/**
	 * Converts the type to the Quantification CHP type.
	 * 
	 * @param chip
	 *          The pointer to the CHP data object.
	 * @return The Quantification CHP data type or NULL if not compatible.
	 */
	public static FusionCHPQuantificationData fromBase(FusionCHPData chip) {
		if (chip == null) {
			return null;
		}
		String chipName = chip.getClass().getName();
		String genName = FusionCHPQuantificationData.class.getName();
		if (chipName.compareTo(genName) == 0) {
			return (FusionCHPQuantificationData)chip;
		}
		return null;
	}
};
