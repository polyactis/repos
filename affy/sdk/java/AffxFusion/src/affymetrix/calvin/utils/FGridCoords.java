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

/** Defines a floating point grid coords */
public class FGridCoords {

	/** Creates a new instance of FGridCoords */
	public FGridCoords() {
		clear();
	}

	/** Clears the members. */
	private void clear() {
		upperleft = new FPoint();
		upperright = new FPoint();
		lowerright = new FPoint();
		lowerleft = new FPoint();
	}

	/**
	 * Cast constructor from FRegion.
	 * 
	 * @param r
	 *          The region of 4 coordinates.
	 */
	public FGridCoords(FRegion r) {
		upperleft = new FPoint(r.get(RectanglePositions.UpperLeft));
		upperright = new FPoint(r.get(RectanglePositions.UpperRight));
		lowerright = new FPoint(r.get(RectanglePositions.LowerRight));
		lowerleft = new FPoint(r.get(RectanglePositions.LowerLeft));
	}

	/**
	 * Return a region
	 * 
	 * @return A region defined by the grid coordinates.
	 */
	public FRegion getRegion() {
		FRegion r = new FRegion();
		r.add(new FPoint(upperleft));
		r.add(new FPoint(upperright));
		r.add(new FPoint(lowerright));
		r.add(new FPoint(lowerleft));
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
	private FPoint upperleft;

	/** Gets the upper left coordinate */
	public FPoint getUpperLeft() {
		return upperleft;
	}

	/** Sets the upper left coordinate */
	public void setUpperLeft(FPoint pt) {
		upperleft = pt;
	}

	/** The upper right coordinate */
	private FPoint upperright;

	/** Gets the upper Right coordinate */
	public FPoint getUpperRight() {
		return upperright;
	}

	/** Sets the upper Right coordinate */
	public void setUpperRight(FPoint pt) {
		upperright = pt;
	}

	/** The lower right coordinate */
	private FPoint lowerright;

	/** Gets the upper left coordinate */
	public FPoint getLowerRight() {
		return lowerright;
	}

	/** Sets the upper left coordinate */
	public void setLowerRight(FPoint pt) {
		lowerright = pt;
	}

	/** The lower left coordinate */
	private FPoint lowerleft;

	/** Gets the lower left coordinate */
	public FPoint getLowerLeft() {
		return lowerleft;
	}

	/** Sets the lower left coordinate */
	public void setLowerLeft(FPoint pt) {
		lowerleft = pt;
	}

}
