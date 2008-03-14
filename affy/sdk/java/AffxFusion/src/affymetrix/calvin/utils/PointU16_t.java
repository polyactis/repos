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

/** Defines an shortegral location */
public class PointU16_t {

	/** Creates a new instance of Poshort */
	public PointU16_t() {
		x = 0;
		y = 0;
	}

	/**
	 * Creates a new instance of Poshort
	 * 
	 * @param pt
	 *          The poshort to copy.
	 */
	public PointU16_t(PointU16_t pt) {
		x = pt.getX();
		y = pt.getY();
	}

	/** The x coordinate */
	private short x;

	/** Gets the x coordinate. */
	public short getX() {
		return x;
	}

	/** Sets the x coordinate. */
	public void setX(short val) {
		x = val;
	}

	/** The y coordinate */
	private short y;

	/** Gets the y coordinate. */
	public short getY() {
		return y;
	}

	/** Sets the y coordinate. */
	public void setY(short val) {
		y = val;
	}

	/** Equality test */
	public boolean equals(PointU16_t lhs) {
		return ((x == lhs.x) && (y == lhs.y));
	}

}
