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

import affymetrix.portability.DataSizes;

/** Defines information about a background zone. */
public class BackgroundZoneType {

	/** The X coordinate of the center of the zone. */
	private float centerX;

	/**
	 * Gets the X coordinate of the center of the zone.
	 * 
	 * @return The X coordinate.
	 */
	public float getCenterX() {
		return centerX;
	}

	/**
	 * Sets the X coordinate of the center of the zone.
	 * 
	 * @param x
	 *          The X coordinate.
	 */
	public void setCenterX(float x) {
		centerX = x;
	}

	/** The Y coordinate of the center of the zone. */
	private float centerY;

	/**
	 * Gets the Y coordinate of the center of the zone.
	 * 
	 * @return The Y coordinate.
	 */
	public float getCenterY() {
		return centerY;
	}

	/**
	 * Sets the Y coordinate of the center of the zone.
	 * 
	 * @param y
	 *          The Y coordinate.
	 */
	public void setCenterY(float y) {
		centerY = y;
	}

	/** The background value for the zone. */
	private float background;

	/**
	 * Gets the background.
	 * 
	 * @return The background value for the zone.
	 */
	public float getBackground() {
		return background;
	}

	/**
	 * Sets the background.
	 * 
	 * @param b
	 *          The background value for the zone.
	 */
	public void setBackground(float b) {
		background = b;
	}

	/** Creates a new instance of BackgroundZoneType */
	public BackgroundZoneType() {
		centerX = 0.0f;
		centerY = 0.0f;
		background = 0.0f;
	}

	/**
	 * Creates a new instance of BackgroundZoneType with values.
	 * 
	 * @param z
	 *          The background zone.
	 */
	public BackgroundZoneType(BackgroundZoneType z) {
		centerX = z.getCenterX();
		centerY = z.getCenterY();
		background = z.getBackground();
	}

	/**
	 * Creates a new instance of BackgroundZoneType with values.
	 * 
	 * @param x
	 *          The x coordinate of the zone.
	 * @param y
	 *          The y coordinate of the zone.
	 * @param b
	 *          The background value.
	 */
	public BackgroundZoneType(float x, float y, float b) {
		centerX = x;
		centerY = y;
		background = b;
	}

	/** The number of bytes used to store the zone information in the CHP file. */
	public static final int ZONE_INFO_TYPE_SIZE = 3 * DataSizes.FLOAT_SIZE;
}
