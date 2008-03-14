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

package affymetrix.gcos.chp;

import java.util.ArrayList;
import java.util.List;

/** Information about all background zones on an array. */
public class BackgroundZoneInfo {

	/** A vector of BackgroundZoneType objects. */
	private List<BackgroundZoneType> zones = new ArrayList<BackgroundZoneType>();

	/**
	 * Gets the number of zones.
	 * 
	 * @return The number of zones on the array.
	 */
	public int getNumberZones() {
		return zones.size();
	}

	/**
	 * Adds a zone
	 * 
	 * @param zone
	 *          The zone to add.
	 */
	public void addZone(BackgroundZoneType zone) {
		zones.add(zone);
	}

	/**
	 * Gets a zone by index.
	 * 
	 * @param index
	 *          The zero based index to the zone of interest.
	 * @return The zone.
	 */
	public BackgroundZoneType getZone(int index) {
		return zones.get(index);
	}

	/** The smoothing factor used to calculate the background values. */
	private float smoothFactor;

	/**
	 * Gets the smoothing factor.
	 * 
	 * @return The smoothing factor.
	 */
	public float getSmoothFactor() {
		return smoothFactor;
	}

	/**
	 * Sets the smoothing factor.
	 * 
	 * @param sf
	 *          The smoothing factor.
	 */
	public void setSmoothFactor(float sf) {
		smoothFactor = sf;
	}

	/** Creates a new instance of BackgroundZoneInfo */
	public BackgroundZoneInfo() {
		smoothFactor = 0.0f;
	}

	/** Clears the members. */
	public void clear() {
		smoothFactor = 0.0f;
		zones.clear();
	}
}
