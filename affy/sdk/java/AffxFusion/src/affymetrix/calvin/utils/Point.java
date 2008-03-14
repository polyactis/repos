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

/** Defines an integral location */
public class Point {

	/** Creates a new instance of Point */
	public Point() {
		x = 0;
		y = 0;
	}

	/**
	 * Creates a new instance of Point
	 * 
	 * @param pt
	 *          The point to copy.
	 */
	public Point(Point pt) {
		x = pt.getX();
		y = pt.getY();
	}

	/** The x coordinate */
	private int x;

	/** Gets the x coordinate. */
	public int getX() {
		return x;
	}

	/** Sets the x coordinate. */
	public void setX(int val) {
		x = val;
	}

	/** The y coordinate */
	private int y;

	/** Gets the y coordinate. */
	public int getY() {
		return y;
	}

	/** Sets the y coordinate. */
	public void setY(int val) {
		y = val;
	}

	/** Equality test */
	public boolean equals(Point lhs) {
		return ((x == lhs.x) && (y == lhs.y));
	}

}
