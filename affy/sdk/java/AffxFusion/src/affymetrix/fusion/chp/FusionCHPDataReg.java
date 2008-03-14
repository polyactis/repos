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
import java.util.List;

import affymetrix.calvin.data.GenericData;
import affymetrix.calvin.parsers.GenericFileReader;
import affymetrix.calvin.utils.AffymetrixGuidType;

/** A class used to self register CHP data classes. */
public abstract class FusionCHPDataReg {

	/** A pointer to the first registered CHP reader. */
	private static FusionCHPDataReg head = null;

	/** A pointer to the next registered CHP reader. */
	private FusionCHPDataReg next;

	/** The file type identifiers associated with the CHP files the reader can parse. */
	private List<AffymetrixGuidType> fileTypeIdentifiers;

	/** Constructor */
	public FusionCHPDataReg() {
		next = head;
		head = this;
	}

	/**
	 * Makes an CHP data object.
	 * 
	 * @return The CHP data object.
	 */
	protected abstract FusionCHPData makeObject();

	/*
	 * Sets the file type ids. @param fileTypeIds The identifiers that the CHP object is compatible with.
	 */
	public void setFileTypeIds(List<AffymetrixGuidType> fileTypeIds) {
		fileTypeIdentifiers = fileTypeIds;
	}

	/*
	 * Read the guid from the file.
	 */
	private static AffymetrixGuidType readGuidFromFile(String fileName) {
		if (!new File(fileName).exists()) {
			return null;
		}
		GenericData data = new GenericData();
		GenericFileReader reader = new GenericFileReader();
		try {
			reader.setFilename(fileName);
			reader.readHeader(data, GenericFileReader.ReadHeaderOption.ReadNoDataGroupHeader);
			AffymetrixGuidType guid = new AffymetrixGuidType();
			guid.setGuid(data.getHeader().getGenericDataHdr().getFileTypeId());
			return guid;
		}
		catch (Throwable t) {
			AffymetrixGuidType guid = new AffymetrixGuidType();
			return guid;
		}
	}

	/**
	 * Reads the contents of a CHP file.
	 * 
	 * @param fileName
	 *          The full path to the input CHP file.
	 * @return A pointer to the CHP data object. NULL if the read failed.
	 */
	public static FusionCHPData read(String fileName) {
		AffymetrixGuidType id = readGuidFromFile(fileName);
		if (id == null) {
			return null;
		}
		FusionCHPData chp = createObject(id);
		if (chp != null) {
			chp.setFileName(fileName);
			if (!chp.read()) {
				chp = null;
			}
		}
		return chp;
	}

	/**
	 * Reads the header of a CHP file.
	 * 
	 * @param fileName
	 *          The full path to the input CHP file.
	 * @return A pointer to the CHP data object. NULL if the read failed.
	 */
	public static FusionCHPData readHeader(String fileName) {
		AffymetrixGuidType id = readGuidFromFile(fileName);
		if (id == null) {
			return null;
		}
		FusionCHPData chp = FusionCHPDataReg.createObject(id);
		if (chp != null) {
			chp.setFileName(fileName);
			if (chp.readHeader() == false) {
				chp = null;
			}
		}
		return chp;
	}

	/**
	 * Creates a CHP reading object.
	 * 
	 * @param fileTypeId
	 *          The file type in the CHP file.
	 * @return The CHP object, NULL if not able to read the file.
	 */
	private static FusionCHPData createObject(AffymetrixGuidType fileTypeId) {
		// Find the matching CHP data object.
		for (FusionCHPDataReg p = head; p != null; p = p.next) {
			List<AffymetrixGuidType> ids = p.fileTypeIdentifiers;
			for (int i = 0; i < ids.size(); i++) {
				AffymetrixGuidType id = ids.get(i);
				if (id.equals(fileTypeId)) {
					FusionCHPData chp = p.makeObject();
					chp.setFileTypeIdentifiers(p.fileTypeIdentifiers);
					chp.setFileTypeIdentifier(fileTypeId);
					if (chp != null) {
						return chp;
					}
				}
			}
		}
		// Go back now and find the generic reader (the one with no file type identifiers)
		for (FusionCHPDataReg p = head; p != null; p = p.next) {
			if (p.fileTypeIdentifiers.size() == 0) {
				FusionCHPData chp = p.makeObject();
				chp.setFileTypeIdentifier(fileTypeId);
				return chp;
			}
		}
		return null;
	}

}
