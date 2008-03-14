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

/** Stores the data for a single genomic data item. */
public class CHPTilingEntry {

	/** The genomic position. */
	private int position;

	/** The genomic position. */
	public int getPosition() {
		return position;
	}

	/** The genomic position. */
	public void setPosition(int p) {
		position = p;
	}

	/** The value associated with the position. */
	private float value;

	/** The value associated with the position. */
	public float getValue() {
		return value;
	}

	/** The value associated with the position. */
	public void setValue(float v) {
		value = v;
	}

	/** Creates a new instance of CHPTilingEntry */
	public CHPTilingEntry() {
		position = 0;
		value = 0.0f;
	}
}
