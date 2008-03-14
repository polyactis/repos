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

/** Defines an integral rectagle */
public class GridCoords {

	/** Creates a new instance of GridCoords */
	public GridCoords() {
		clear();
	}

	/** Clears the members. */
	private void clear() {
		upperleft = new Point();
		upperright = new Point();
		lowerright = new Point();
		lowerleft = new Point();
	}

	/**
	 * Cast constructor from FRegion.
	 * 
	 * @param r
	 *          The region of 4 coordinates.
	 */
	public GridCoords(Region r) {
		upperleft = new Point(r.get(RectanglePositions.UpperLeft));
		upperright = new Point(r.get(RectanglePositions.UpperRight));
		lowerright = new Point(r.get(RectanglePositions.LowerRight));
		lowerleft = new Point(r.get(RectanglePositions.LowerLeft));
	}

	/**
	 * Return a region
	 * 
	 * @return A region defined by the grid coordinates.
	 */
	public Region getRegion() {
		Region r = new Region();
		r.add(new Point(upperleft));
		r.add(new Point(upperright));
		r.add(new Point(lowerright));
		r.add(new Point(lowerleft));
		return r;
	}

	/** Tests if the rectangle is empty */
	public boolean isEmpty() {
		if ((upperleft == null) || (upperright == null) || (lowerright == null) || (lowerleft == null)) {
			return true;
		}
		else if (upperleft.equals(upperright) && lowerleft.equals(lowerright) && upperleft.equals(lowerleft)) {
			return true;
		}
		return false;
	}

	/** The upper left coordinate */
	private Point upperleft;

	/** Gets the upper left coordinate */
	public Point getUpperLeft() {
		return upperleft;
	}

	/** Sets the upper left coordinate */
	public void setUpperLeft(Point pt) {
		upperleft = pt;
	}

	/** The upper right coordinate */
	private Point upperright;

	/** Gets the upper Right coordinate */
	public Point getUpperRight() {
		return upperright;
	}

	/** Sets the upper Right coordinate */
	public void setUpperRight(Point pt) {
		upperright = pt;
	}

	/** The lower right coordinate */
	private Point lowerright;

	/** Gets the upper left coordinate */
	public Point getLowerRight() {
		return lowerright;
	}

	/** Sets the upper left coordinate */
	public void setLowerRight(Point pt) {
		lowerright = pt;
	}

	/** The lower left coordinate */
	private Point lowerleft;

	/** Gets the lower left coordinate */
	public Point getLowerLeft() {
		return lowerleft;
	}

	/** Sets the lower left coordinate */
	public void setLowerLeft(Point pt) {
		lowerleft = pt;
	}
}
