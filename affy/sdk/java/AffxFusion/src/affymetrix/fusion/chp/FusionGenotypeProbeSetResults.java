//////////////////////////////////////////////////////////////////
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

package affymetrix.fusion.chp;

import affymetrix.calvin.data.CHPGenotypeEntry;
import affymetrix.gcos.chp.GenotypeProbeSetResults;

/** Stores results for a genotyping probe set. */
public class FusionGenotypeProbeSetResults {

	/** The AA allele call */
	public static final byte ALLELE_A_CALL = 6;

	/** The BB allele call */
	public static final byte ALLELE_B_CALL = 7;

	/** The AB allele call */
	public static final byte ALLELE_AB_CALL = 8;

	/** The no call allele call */
	public static final byte ALLELE_NO_CALL = 11;

	/** The GCOS probe set results object. */
	private GenotypeProbeSetResults gcosResult;

	/**
	 * Sets the GCOS probe set results object.
	 * 
	 * @param r
	 *          The GCOS probe set results object.
	 */
	public void setGCOSObject(GenotypeProbeSetResults r) {
		gcosResult = r;
	}

	/** The calvin probe set results object. */
	private CHPGenotypeEntry calvinResult;

	/**
	 * Sets the calvin object.
	 * 
	 * @param r
	 *          The calvin probe set object.
	 */
	public void setCalvinObject(CHPGenotypeEntry r) {
		calvinResult = r;
	}

	/** Clears the members. */
	public void clear() {
		gcosResult = null;
		calvinResult = null;
	}

	/**
	 * Gets the allele call
	 * 
	 * @return The allele call.
	 */
	public byte getAlleleCall() {
		if (gcosResult != null) {
			return gcosResult.getAlleleCall();
		}
		else if (calvinResult != null) {
			return calvinResult.getCall();
		}
		return 0;
	}

	/**
	 * Gets the confidence
	 * 
	 * @return The confidence.
	 */
	public float getConfidence() {
		if (gcosResult != null) {
			return gcosResult.getConfidence();
		}
		else if (calvinResult != null) {
			return calvinResult.getConfidence();
		}
		return 0.0f;
	}

	/**
	 * Gets the RAS1 value
	 * 
	 * @return The RAS1 value.
	 */
	public float getRAS1() {
		if (gcosResult != null) {
			return gcosResult.getRAS1();
		}
		else if (calvinResult != null) {
			return calvinResult.getRAS1();
		}
		return 0.0f;
	}

	/**
	 * Gets the RAS2 value
	 * 
	 * @return The RAS2 value.
	 */
	public float getRAS2() {
		if (gcosResult != null) {
			return gcosResult.getRAS2();
		}
		else if (calvinResult != null) {
			return calvinResult.getRAS2();
		}
		return 0.0f;
	}

	/**
	 * Gets the AA p-value
	 * 
	 * @return The AA p-value.
	 */
	public float getPValue_AA() {
		if (gcosResult != null) {
			return gcosResult.getPValue_AA();
		}
		else if (calvinResult != null) {
			return calvinResult.getAACall();
		}
		return 0.0f;
	}

	/**
	 * Gets the AB p-value
	 * 
	 * @return The AB p-value.
	 */
	public float getPValue_AB() {
		if (gcosResult != null) {
			return gcosResult.getPValue_AB();
		}
		else if (calvinResult != null) {
			return calvinResult.getABCall();
		}
		return 0.0f;
	}

	/**
	 * Gets the BB p-value
	 * 
	 * @return The BB p-value.
	 */
	public float getPValue_BB() {
		if (gcosResult != null) {
			return gcosResult.getPValue_BB();
		}
		else if (calvinResult != null) {
			return calvinResult.getBBCall();
		}
		return 0.0f;
	}

	/**
	 * Gets the NoCall p-value
	 * 
	 * @return The NoCall p-value.
	 */
	public float getPValue_NoCall() {
		if (gcosResult != null) {
			return gcosResult.getPValue_NoCall();
		}
		else if (calvinResult != null) {
			return calvinResult.getNoCall();
		}
		return 0.0f;
	}

	/**
	 * Returns a string representation of the allele call.
	 * 
	 * @return The allele call as a string.
	 */
	public String getAlleleCallString() {
		switch (getAlleleCall()) {
		case ALLELE_A_CALL:
			return "A";

		case ALLELE_B_CALL:
			return "B";

		case ALLELE_AB_CALL:
			return "AB";

		default:
			return "No Call";
		}
	}

	/** Creates a new instance of GenotypeProbeSetResults */
	public FusionGenotypeProbeSetResults() {
		clear();
	}

}
