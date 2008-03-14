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

/** This class stores a zone's background value */
public class CHPUniversalEntry {

	/** Creates a new instance of CHPUniversalEntry */
	public CHPUniversalEntry() {
		clear();
	}

	/** Creates a new instance of CHPUniversalEntry */
	public CHPUniversalEntry(float bg) {
		background = bg;
	}

	/** Creates a new instance of CHPUniversalEntry */
	public CHPUniversalEntry(CHPUniversalEntry entry) {
		background = entry.getBackground();
	}

	/** Clears the members. */
	public void clear() {
		background = 0.0f;
	}

	/** The background value. */
	private float background;

	public float getBackground() {
		return background;
	}

	public void setBackground(float p) {
		background = p;
	}
}
