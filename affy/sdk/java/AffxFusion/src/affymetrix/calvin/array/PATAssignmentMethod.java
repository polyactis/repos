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

/** The method the probe array type was assigned. */
public class PATAssignmentMethod {

	/** No probe array type assignment. */
	public static final String PAT_ASSIGNMENT_NONE = "None";

	/** Affy barcode selected probe array type. */
	public static final String PAT_ASSIGNMENT_BARCODE = "AffyBarcode";

	/** User selected probe array type. */
	public static final String PAT_ASSIGNMENT_USER_SELECTED = "UserSelected";

	/** Other assignment. */
	public static final String PAT_ASSIGNMENT_OTHER = "Other";

	/** The method the probe array type was assigned. */
	public enum PATAssignmentMethodType {
		NoAssignment, // Unknown
		AffyBarcodeAssignment, // User entered an Affy barcode.
		UserSelectedAssignment, // User selected.
		OtherAssignment
		// Other method.
	}

	private PATAssignmentMethod() {
	}

	/**
	 * Converts a string to the probe array assignment type.
	 * 
	 * @param str
	 *          The string representation.
	 */
	public static PATAssignmentMethodType toPATAssignmentMethodType(String str) {
		if (str.equals(PAT_ASSIGNMENT_BARCODE)) {
			return PATAssignmentMethodType.AffyBarcodeAssignment;
		}
		else if (str.equals(PAT_ASSIGNMENT_USER_SELECTED)) {
			return PATAssignmentMethodType.UserSelectedAssignment;
		}
		else if (str.equals(PAT_ASSIGNMENT_OTHER)) {
			return PATAssignmentMethodType.OtherAssignment;
		}
		else {
			return PATAssignmentMethodType.NoAssignment;
		}
	}

	/**
	 * Converts the probe array assignment type to a string.
	 * 
	 * @return The string representation.
	 */
	public static String toString(PATAssignmentMethodType t) {
		switch (t) {
		case NoAssignment:
			return PAT_ASSIGNMENT_NONE;

		case AffyBarcodeAssignment:
			return PAT_ASSIGNMENT_BARCODE;

		case UserSelectedAssignment:
			return PAT_ASSIGNMENT_USER_SELECTED;

		case OtherAssignment:
			return PAT_ASSIGNMENT_OTHER;
		}
		return null;
	}

}
