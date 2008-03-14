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

package affymetrix.calvin.array;

/** Defines constants for the media an array is housed. */
public class ArrayMedia {

	/** The GeneChip cartridge. */
	public static final String MEDIA_TYPE_CARTRIDGE = "Cartridge";

	/** A 96 well plate or strip. */
	public static final String MEDIA_TYPE_PLATE_OR_STRIP = "PlateOrStrip";
	
	/*! The type of media of the array. */
	public enum ArrayMediaType
	{
		CartridgeMedia,		/*! A GeneChip cartridge. */
		PlateOrStripMedia	/*! A 96 well plate or peg strip. */
	}
	
	/** Ctor */
	private ArrayMedia() {
	}

	/**
	 * Converts the media type to a string.
	 * 
	 * @return The string representation.
	 */
	public static String toString(ArrayMediaType t) {
		switch (t) {
		case CartridgeMedia:
			return MEDIA_TYPE_CARTRIDGE;

		case PlateOrStripMedia:
			return MEDIA_TYPE_PLATE_OR_STRIP;
		}
		return null;
	}

	/**
	 * Creates a new instance given a string representation of the media.
	 * 
	 * @param str
	 *          The string representation.
	 */
	public static ArrayMediaType toArrayMediaType(String str) {
		if (str.equals(MEDIA_TYPE_PLATE_OR_STRIP)) {
			return ArrayMediaType.PlateOrStripMedia;
		}
		else {
			return ArrayMediaType.CartridgeMedia;
		}
	}
}
