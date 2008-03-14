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

package affymetrix.calvin.data;

/** This class stores a zone's background value */
public class CHPBackgroundZone {

	/** Creates a new instance of CHPBackgroundZone */
	public CHPBackgroundZone() {
		clear();
	}

	/** Creates a new instance of CHPBackgroundZone */
	public CHPBackgroundZone(float x, float y, float bg, float smooth) {
		centerX = x;
		centerY = y;
		background = bg;
		smoothFactor = smooth;
	}

	/** Creates a new instance of CHPBackgroundZone */
	public CHPBackgroundZone(CHPBackgroundZone zn) {
		centerX = zn.getCenterX();
		centerY = zn.getCenterY();
		background = zn.getBackground();
		smoothFactor = zn.getSmoothFactor();
	}

	/** The X coordinate of the center of the zone. */
	private float centerX;

	/** The Y coordinate of the center of the zone. */
	private float centerY;

	/** The zone's background value */
	private float background;

	/** The smoothing factor used to calculate the zone backgrounds */
	private float smoothFactor;

	public void clear() {
		centerX = 0.0f;
		centerY = 0.0f;
		background = 0.0f;
		smoothFactor = 0.0f;
	}

	public float getCenterX() {
		return centerX;
	}

	public float getCenterY() {
		return centerY;
	}

	public float getBackground() {
		return background;
	}

	public float getSmoothFactor() {
		return smoothFactor;
	}

	public void setCenterX(float p) {
		centerX = p;
	}

	public void setCenterY(float p) {
		centerY = p;
	}

	public void setBackground(float p) {
		background = p;
	}

	public void setSmoothFactor(float p) {
		smoothFactor = p;
	}
}
