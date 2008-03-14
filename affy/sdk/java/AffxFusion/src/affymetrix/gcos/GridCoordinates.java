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

/** Stores the coordinates for a global grid (4 corner grid). */
public class GridCoordinates {

	/** The upper left corner. */
	private CoordinatePoint upperleft;

	/**
	 * Gets the upper left corner.
	 * 
	 * @return The coordinate.
	 */
	public CoordinatePoint getUpperLeft() {
		return upperleft;
	}

	/**
	 * Sets the upper left corner.
	 * 
	 * @param coord
	 *          The coordinate.
	 */
	public void setUpperLeft(CoordinatePoint coord) {
		upperleft = new CoordinatePoint(coord);
	}

	/** The upper right corner. */
	private CoordinatePoint upperright;

	/**
	 * Gets the upper right corner.
	 * 
	 * @return The coordinate.
	 */
	public CoordinatePoint getUpperRight() {
		return upperright;
	}

	/**
	 * Sets the upper right corner.
	 * 
	 * @param coord
	 *          The coordinate.
	 */
	public void setUpperRight(CoordinatePoint coord) {
		upperright = new CoordinatePoint(coord);
	}

	/** The lower right corner. */
	private CoordinatePoint lowerright;

	/**
	 * Gets the lower right corner.
	 * 
	 * @return The coordinate.
	 */
	public CoordinatePoint getLowerRight() {
		return lowerright;
	}

	/**
	 * Sets the lower right corner.
	 * 
	 * @param coord
	 *          The coordinate.
	 */
	public void setLowerRight(CoordinatePoint coord) {
		lowerright = new CoordinatePoint(coord);
	}

	/** The lower left corner. */
	private CoordinatePoint lowerleft;

	/**
	 * Gets the lower left corner.
	 * 
	 * @return The coordinate.
	 */
	public CoordinatePoint getLowerLeft() {
		return lowerleft;
	}

	/**
	 * Sets the lower left corner.
	 * 
	 * @param coord
	 *          The coordinate.
	 */
	public void setLowerLeft(CoordinatePoint coord) {
		lowerleft = new CoordinatePoint(coord);
	}

	/** Creates a new instance of GridCoordinates */
	public GridCoordinates() {
		upperleft = null;
		upperright = null;
		lowerleft = null;
		lowerright = null;
	}

	/** Creates a new instance of GridCoordinates */
	public GridCoordinates(GridCoordinates grid) {
		upperleft = new CoordinatePoint(grid.getUpperLeft());
		upperright = new CoordinatePoint(grid.getUpperRight());
		lowerleft = new CoordinatePoint(grid.getLowerLeft());
		lowerright = new CoordinatePoint(grid.getLowerRight());
	}

}
