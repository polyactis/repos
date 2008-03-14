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

package affymetrix.calvin.data;

import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.portability.UInt;

/** This class defines a data container for the generic file header */
public class FileHeader {

	/** The magic number for calvin files. */
	static public final byte MAGIC_NUM = 59;

	/** The version number for calvin files. */
	static public final byte VERSION = 1;

	/** filename */
	private String filename;

	private byte magic;

	private byte version;

	private List<DataGroupHeader> dataGroupHdrs = new ArrayList<DataGroupHeader>();

	private GenericDataHeader genericHdr = new GenericDataHeader();

	/** Number of data dataGroups in the file */
	private UInt numDataGroups = new UInt();

	/** Position of the first DataGroup. */
	private UInt firstDataGroupFilePos = UInt.ZERO;

	/** Creates a new instance of FileHeader */
	public FileHeader() {
		magic = MAGIC_NUM;
		version = VERSION;
	}

	public void clear() {
		dataGroupHdrs.clear();
		genericHdr.clear();
	}

	public void setFilename(String p) {
		filename = p;
	}

	public String getFilename() {
		return filename;
	}

	public byte getMagicNumber() {
		return magic;
	}

	public byte getVersion() {
		return version;
	}

	/** Get the number of DataGroupHeaders added */
	public int getDataGroupCnt() {
		return dataGroupHdrs.size();
	}

	public void addDataGroupHdr(DataGroupHeader p) {
		dataGroupHdrs.add(p);
	}

	/** Get a DataGroupHeader by index. Max index < GetDataGroupCnt */
	public DataGroupHeader getDataGroup(int index) {
		if (dataGroupHdrs.size() > 0) {
			return dataGroupHdrs.get(index);
		}
		return null;
	}

	public void setGenericDataHdr(GenericDataHeader p) {
		genericHdr = p;
	}

	public GenericDataHeader getGenericDataHdr() {
		return genericHdr;
	}

	/**
	 * Finds a DataGroupHeader by name.
	 * 
	 * @param name
	 *          The name of the DataGroup
	 * @return A pointer to the DataGroupHeader. If not found, the return is 0.
	 */
	public DataGroupHeader findDataGroupHeader(String name) {
		int n = getDataGroupCnt();
		for (int i = 0; i < n; i++) {
			DataGroupHeader h = getDataGroup(i);
			if (name.equals(h.getName())) {
				return h;
			}
		}
		return null;
	}

	/** Get the number of DataGroups in a file. */
	public UInt getNumDataGroups() {
		return numDataGroups;
	}

	/** Set the number of DataGroups. Set when reading a file */
	public void setNumDataGroups(UInt value) throws UnsignedOutOfLimitsException {
		numDataGroups.set(value.toLong());
	}

	/** Get the file position to the first DataGroup */
	public UInt getFirstDataGroupFilePos() {
		return firstDataGroupFilePos;
	}

	/** Set the file postion to the first DataGroup. Method should be protected. It is set when the object is being read. */
	public void setFirstDataGroupFilePos(UInt value) {
		firstDataGroupFilePos = value;
	}
}
