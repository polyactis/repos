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

import java.util.ArrayList;
import java.util.List;

/** Defines a region with integral point coordinates */
public class Region {

	/** A vector of points */
	private List<Point> pts = new ArrayList<Point>();

	/** Creates a new instance of Region */
	public Region() {
		clear();
	}

	/** Clears the region. */
	public void clear() {
		pts.clear();
	}

	/** Gets the size of the Point vector. */
	public int size() {
		return pts.size();
	}

	/**
	 * Gets the element at the given index position.
	 * 
	 * @param index
	 *          The zero based index to the vector.
	 * @return The integral point coordinate.
	 */
	public Point get(int index) {
		return pts.get(index);
	}

	/**
	 * Adds a new value to the vector.
	 * 
	 * @param pt
	 *          The new value to add.
	 */
	public void add(Point pt) {
		pts.add(pt);
	}

	/** Equality test */
	public boolean equals(Region lhs) {
		if (lhs.pts.size() == pts.size()) {
			int sz = pts.size();
			for (int i = 0; i < sz; i++) {
				if (lhs.get(i).equals(pts.get(i))) {
					continue;
				}
				return false;
			}
			return true;
		}
		return false;
	}
}
