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

package affymetrix.calvin.utils;

import java.io.UnsupportedEncodingException;
import java.rmi.dgc.VMID;
import java.util.Arrays;

/** Provides a type for the Affy GUID. */
public class AffymetrixGuidType {

	public static final int GUID_LENGTH = 54;

	/** The guid. */
	private byte[] guid;

	/** Creates a new instance of AffymetrixGuidType */
	public AffymetrixGuidType() {
		guid = new byte[GUID_LENGTH];
		clear();
	}
	
	public AffymetrixGuidType(boolean generateGuid) {
		this();
		if(generateGuid){
			generateGuid();
		}
	}

	public AffymetrixGuidType(String g) {
		this();
		setGuid(g);
	}

	public AffymetrixGuidType(byte[] g) {
		this();
		setGuid(g);
	}

	/** Gets the guid. */
	public byte[] getGuid() {
		return guid;
	}

	/**
	 * Sets the guid.
	 * 
	 * @param g
	 *          The guid.
	 */
	public void setGuid(byte[] g) {
		if (g.length < GUID_LENGTH) {
			clear();
			System.arraycopy(g, 0, guid, 0, g.length);
		}
		else {
			System.arraycopy(g, 0, guid, 0, GUID_LENGTH);
		}
	}

	public void setGuid(String g) {
		try {
			setGuid(g.getBytes("US-ASCII"));
		}
		catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		}
	}

	public void setGuid(AffymetrixGuidType g) {
		guid = g.getGuid();
	}

	public void clear() {
		for (int i = 0; i < GUID_LENGTH; i++) {
			guid[i] = 0;
		}
	}

	public boolean isEmpty() {
		for (int i = 0; i < GUID_LENGTH; i++) {
			if (guid[i] != 0) {
				return false;
			}
		}
		return true;
	}

	public boolean equals(AffymetrixGuidType g) {
		return Arrays.equals(guid, g.getGuid());
	}

	public void generateGuid() {
		String s = new VMID().toString();
		s = s.replace(':', '-');
		setGuid(s);
	}

	@Override
	public String toString() {
		try {
			return new String(guid, IOUtils.ASCII_CHARSET);
		}
		catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		}
		return null;
	}
}
