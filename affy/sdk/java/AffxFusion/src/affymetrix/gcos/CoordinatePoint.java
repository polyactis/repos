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

package affymetrix.gcos;

/** Defines a class to hold a coordinate. */
public class CoordinatePoint {

	/** The X coordinate. */
	private int x;

	/**
	 * Gets the X coordinate.
	 * 
	 * @return The X coordinate.
	 */
	public int getX() {
		return x;
	}

	/**
	 * Sets the X coordinate.
	 * 
	 * @param xcoord
	 *          The X coordinate.
	 */
	public void setX(int xcoord) {
		x = xcoord;
	}

	/** The Y coordinate. */
	private int y;

	/**
	 * Gets the Y coordinate.
	 * 
	 * @return The Y coordinate.
	 */
	public int getY() {
		return y;
	}

	/**
	 * Sets the Y coordinate.
	 * 
	 * @param ycoord
	 *          The Y coordinate.
	 */
	public void setY(int ycoord) {
		y = ycoord;
	}

	/** Creates a new instance of CoordinatePoint */
	public CoordinatePoint() {
		x = 0;
		y = 0;
	}

	/** Creates a new instance of CoordinatePoint */
	public CoordinatePoint(CoordinatePoint coord) {
		x = coord.getX();
		y = coord.getY();
	}

}
