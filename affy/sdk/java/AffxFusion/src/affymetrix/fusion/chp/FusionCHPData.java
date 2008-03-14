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

import java.util.List;

import affymetrix.calvin.utils.AffymetrixGuidType;

/** A base class for all CHP data objects. */
public abstract class FusionCHPData {

	/** The CHP file name. */
	protected String filename;

	/**
	 * Sets the file name.
	 * 
	 * @param str
	 *          The file name.
	 */
	protected void setFileName(String str) {
		filename = str;
	}

	/**
	 * Gets the file name.
	 * 
	 * @return The file name.
	 */
	public String getFileName() {
		return filename;
	}

	/**
	 * Reads the contents of the file.
	 * 
	 * @return True if successfully read.
	 */
	protected abstract boolean read();

	/**
	 * Reads the header.
	 * 
	 * @return True if successfully read.
	 */
	protected abstract boolean readHeader();

	/**
	 * Get the id of the file (only valid for Command Console "calvin" files)
	 * 
	 * @return The unique file id.
	 */
	public abstract AffymetrixGuidType getFileId();

	/** The file type identifiers associated with the CHP files the reader can parse. */
	protected List<AffymetrixGuidType> fileTypeIdentifiers;

	/** The actual file type identifier in the file. */
	protected AffymetrixGuidType fileTypeIdentifier;

	/**
	 * Gets the file type identifiers associated with the CHP files the reader can parse.
	 * 
	 * @return The ids
	 */
	public List<AffymetrixGuidType> getFileTypeIdentifiers() {
		return fileTypeIdentifiers;
	}

	/**
	 * Sets the file type identifiers associated with the CHP files the reader can parse.
	 * 
	 * @param ids
	 *          The ids
	 */
	public void setFileTypeIdentifiers(List<AffymetrixGuidType> ids) {
		fileTypeIdentifiers = ids;
	}

	/**
	 * Gets the file type identifier in the file (blank for GCOS files).
	 * 
	 * @return The id of the file.
	 */
	public AffymetrixGuidType getFileTypeIdentifier() {
		return fileTypeIdentifier;
	}

	/**
	 * Sets the file type identifier in the file (blank for GCOS files).
	 * 
	 * @param id
	 *          The id of the file.
	 */
	public void setFileTypeIdentifier(AffymetrixGuidType id) {
		fileTypeIdentifier = id;
	}
}
