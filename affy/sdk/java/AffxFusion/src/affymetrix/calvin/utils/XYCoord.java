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

/** This is a class for holding x-y coordinates. */
public class XYCoord {

	/** Constructor */
	public XYCoord() {
		xCoord = 0;
		yCoord = 0;
	}

	/**
	 * Constructor
	 * 
	 * @param x
	 *          The x coordinate.
	 * @param y
	 *          The y coordinate.
	 */
	public XYCoord(short x, short y) {
		xCoord = x;
		yCoord = y;
	}

	/**
	 * Constructor
	 * 
	 * @param coord
	 *          The coordinate.
	 */
	public XYCoord(XYCoord coord) {
		xCoord = coord.getX();
		yCoord = coord.getY();
	}

	/** x-coordinate */
	private short xCoord;

	/** Gets the X coordinate. */
	public short getX() {
		return xCoord;
	}

	/** Gets the X coordinate. */
	public void setX(short x) {
		xCoord = x;
	}

	/** y-coordinate */
	private short yCoord;

	/** Gets the Y coordinate. */
	public short getY() {
		return yCoord;
	}

	/** Gets the Y coordinate. */
	public void setY(short y) {
		yCoord = y;
	}

	/** equality operator */
	public boolean equals(XYCoord p) {
		return ((xCoord == p.getX()) && (yCoord == p.getY()));
	}

	/** less than operator */
	public boolean lessThan(XYCoord p) {
		return (yCoord < p.getY() ? true : (((yCoord == p.getY()) && (xCoord < p.getX())) ? true : false));
	}
}
