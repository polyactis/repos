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

/** Defines a floating-point location */
public class FPoint {

	/** Creates a new instance of FPoint */
	public FPoint() {
		x = 0.0f;
		y = 0.0f;
	}

	/**
	 * Creates a new instance of FPoint
	 * 
	 * @param pt
	 *          The point to copy.
	 */
	public FPoint(FPoint pt) {
		x = pt.getX();
		y = pt.getY();
	}

	/** The x coordinate */
	private float x;

	/** Gets the x coordinate. */
	public float getX() {
		return x;
	}

	/** Sets the x coordinate. */
	public void setX(float val) {
		x = val;
	}

	/** The y coordinate */
	private float y;

	/** Gets the y coordinate. */
	public float getY() {
		return y;
	}

	/** Sets the y coordinate. */
	public void setY(float val) {
		y = val;
	}

	/** Equality test */
	public boolean equals(FPoint lhs) {
		return ((x == lhs.x) && (y == lhs.y));
	}
}
