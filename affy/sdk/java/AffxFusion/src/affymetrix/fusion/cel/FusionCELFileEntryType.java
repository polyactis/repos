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

package affymetrix.fusion.cel;

/** Stores the CEL data for a single feature. */
public class FusionCELFileEntryType {

	/** The intensity. */
	private float intensity;

	/**
	 * Gets the intensity.
	 * 
	 * @return The intensity.
	 */
	public float getIntensity() {
		return intensity;
	}

	/**
	 * Sets the intensity
	 * 
	 * @param value
	 *          The intensity
	 */
	public void setIntensity(float value) {
		intensity = value;
	}

	/** The stdv. */
	private float stdv;

	/**
	 * Gets the stdv.
	 * 
	 * @return The stdv.
	 */
	public float getStdv() {
		return stdv;
	}

	/**
	 * Sets the stdv
	 * 
	 * @param value
	 *          The stdv
	 */
	public void setStdv(float value) {
		stdv = value;
	}

	/** The pixels. */
	private short pixels;

	/**
	 * Gets the pixels.
	 * 
	 * @return The pixels.
	 */
	public short getPixels() {
		return pixels;
	}

	/**
	 * Sets the pixels
	 * 
	 * @param value
	 *          The pixels
	 */
	public void setPixels(short value) {
		pixels = value;
	}

	/** Creates a new instance of CELFileEntryType */
	public FusionCELFileEntryType() {
		intensity = 0.0f;
		stdv = 0.0f;
		pixels = 0;
	}

}
