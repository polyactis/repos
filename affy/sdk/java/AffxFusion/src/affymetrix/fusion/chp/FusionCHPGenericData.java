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

import affymetrix.calvin.data.GenericData;
import affymetrix.calvin.parsers.GenericFileReader;
import affymetrix.calvin.utils.AffymetrixGuidType;

/** This class provides storage and reading capabilities for CHP files */
public class FusionCHPGenericData extends FusionCHPData {

	/** Deallocates any memory used by the class object */
	public void clear() {
		genericData = null;
	}

	/**
	 * Gets the generic data object.
	 * 
	 * @return The generic data object.
	 */
	public GenericData getData() {
		return genericData;
	}

	/** Constructor */
	private FusionCHPGenericData() {
		clear();
	}

	/** The generic file data object. */
	private GenericData genericData;

	/**
	 * Reads the CHP file.
	 * 
	 * @return True if successful.
	 */
	@Override
	protected boolean read() {
		genericData = new GenericData();
		GenericFileReader reader = new GenericFileReader();
		try {
			reader.setFilename(filename);
			reader.open(genericData, GenericFileReader.OpenHint.All);
			return true;
		}
		catch (Throwable t) {
			clear();
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
		return genericData.getFileIdentifier();
	}

	/** A class to register a generic CHP reader. */
	private static class Reg extends FusionCHPDataReg {
		/** Constructor - register the file type. */
		public Reg() {
			super();
			List<AffymetrixGuidType> ids = new ArrayList<AffymetrixGuidType>();
			setFileTypeIds(ids);
		}

		/**
		 * Creates a generic CHP object.
		 * 
		 * @return The generic CHP object.
		 */
		@Override
		public FusionCHPData makeObject() {
			return new FusionCHPGenericData();
		}
	};

	/** The one and only registration object. This registers the class as a CHP reader. */
	private static Reg reg = new Reg();

	/** Register the reader with the system. */
	public static void registerReader() {
		if (FusionCHPGenericData.reg == null) {
			FusionCHPGenericData.reg = new Reg();
		}
	}

	/**
	 * Converts the type to the generic CHP type.
	 * 
	 * @param chip
	 *          The pointer to the CHP data object.
	 * @return The generic CHP data type or NULL if not compatible.
	 */
	public static FusionCHPGenericData fromBase(FusionCHPData chip) {
		if (chip == null) {
			return null;
		}
		String chipName = chip.getClass().getName();
		String genName = FusionCHPGenericData.class.getName();
		if (chipName.compareTo(genName) == 0) {
			return (FusionCHPGenericData)chip;
		}
		return null;
	}

}
