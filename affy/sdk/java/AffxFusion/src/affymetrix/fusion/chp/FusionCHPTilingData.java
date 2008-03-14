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

import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.data.CHPTilingData;
import affymetrix.calvin.data.CHPTilingEntry;
import affymetrix.calvin.data.TilingSequenceData;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.parsers.CHPTilingFileReader;
import affymetrix.calvin.utils.AffymetrixGuidType;

/** This class provides storage and reading capabilities for tiling CHP files */
public class FusionCHPTilingData extends FusionCHPData {

	/** Creates a new instance of FusionCHPTilingData */
	private FusionCHPTilingData() {
		chpData = null;
	}

	/** The CHP object. */
	private CHPTilingData chpData;

	/**
	 * Gets the number of sequence stored in the file.
	 * 
	 * @return The number of sequence stored in the file.
	 */
	public int getNumberSequences() {
		return chpData.getNumberSequences();
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

	/**
	 * Gets the algorithm parameters
	 * 
	 * @return The algoirhtm parameters.
	 */
	public List<ParameterNameValue> getAlgParams() {
		return chpData.getAlgParams();
	}

	/**
	 * Gets the sequence data.
	 * 
	 * @return The data associated with the sequence.
	 */
	public TilingSequenceData getTilingSequenceData() {
		return chpData.getTilingSequenceData();
	}

	/**
	 * Gets the number of entries in a tiling sequence.
	 * 
	 * @param index
	 *          The sequence index.
	 * @return The number of entries in the sequence.
	 */
	public int getTilingSequenceEntryCount(int index) {
		return chpData.getTilingSequenceEntryCount(index);
	}

	/**
	 * Opens a group for reading.
	 * 
	 * @param index
	 *          The index to the sequence.
	 */
	public void openTilingSequenceDataSet(int index) {
		chpData.openTilingSequenceDataSet(index);
	}

	/**
	 * Returns the entry for the given row. The data set must be open.
	 * 
	 * @param row
	 *          The row index.
	 * @return The entry.
	 */
	public CHPTilingEntry getTilingSequenceEntry(int row) throws UnsignedOutOfLimitsException {
		return chpData.getTilingSequenceEntry(row);
	}

	/**
	 * Reads the CHP file.
	 * 
	 * @return True if successful.
	 */
	@Override
	protected boolean read() {
		CHPTilingFileReader reader = new CHPTilingFileReader();
		chpData = new CHPTilingData();
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

	/**
	 * Get the id of the file (only valid for Command Console "calvin" files)
	 * 
	 * @return The unique file id.
	 */
	@Override
	public AffymetrixGuidType getFileId() {
		return chpData.getGenericData().getFileIdentifier();
	}

	/** A class to register the tiling CHP reader. */
	private static class Reg extends FusionCHPDataReg {
		/** Constructor - register the tiling file type. */
		public Reg() {
			super();
			List<AffymetrixGuidType> ids = new ArrayList<AffymetrixGuidType>();
			AffymetrixGuidType guid = new AffymetrixGuidType();
			guid.setGuid(CHPTilingData.CHP_TILING_TYPE);
			ids.add(guid);
			setFileTypeIds(ids);
		}

		/**
		 * Creates a tiling CHP object.
		 * 
		 * @return The tiling CHP object.
		 */
		@Override
		public FusionCHPData makeObject() {
			return new FusionCHPTilingData();
		}
	};

	/** The one and only registration object. This registers the class as a CHP reader. */
	private static Reg reg = null;

	/** Register the reader with the system. */
	public static void registerReader() {
		if (FusionCHPTilingData.reg == null) {
			FusionCHPTilingData.reg = new Reg();
		}
	}

	/**
	 * Converts the type to the tiling CHP type.
	 * 
	 * @param chip
	 *          The pointer to the CHP data object.
	 * @return The tiling CHP data type or NULL if not compatible.
	 */
	public static FusionCHPTilingData fromBase(FusionCHPData chip) {
		if (chip == null) {
			return null;
		}
		String chipName = chip.getClass().getName();
		String genName = FusionCHPTilingData.class.getName();
		if (chipName.compareTo(genName) == 0) {
			return (FusionCHPTilingData)chip;
		}
		return null;
	}
}
