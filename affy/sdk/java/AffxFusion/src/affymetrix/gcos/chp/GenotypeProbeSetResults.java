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

package affymetrix.gcos.chp;

/** Stores results for a genotyping probe set. */
public class GenotypeProbeSetResults {

	/** The AA allele call */
	public static final byte ALLELE_A_CALL = 6;

	/** The BB allele call */
	public static final byte ALLELE_B_CALL = 7;

	/** The AB allele call */
	public static final byte ALLELE_AB_CALL = 8;

	/** The no call allele call */
	public static final byte ALLELE_NO_CALL = 11;

	/** The allele call */
	private byte alleleCall;

	/**
	 * Gets the allele call
	 * 
	 * @return The allele call.
	 */
	public byte getAlleleCall() {
		return alleleCall;
	}

	/**
	 * Sets the allele call
	 * 
	 * @param val
	 *          The allele call.
	 */
	public void setAlleleCall(byte val) {
		alleleCall = val;
	}

	/** The confidence associated with the allele call */
	private float confidence;

	/**
	 * Gets the confidence
	 * 
	 * @return The confidence.
	 */
	public float getConfidence() {
		return confidence;
	}

	/**
	 * Sets the confidence
	 * 
	 * @param val
	 *          The confidence.
	 */
	public void setConfidence(float val) {
		confidence = val;
	}

	/** The relative allele strength for the first block */
	private float ras1;

	/**
	 * Gets the RAS1 value
	 * 
	 * @return The RAS1 value.
	 */
	public float getRAS1() {
		return ras1;
	}

	/**
	 * Sets the RAS1 value
	 * 
	 * @param val
	 *          The RAS1 value.
	 */
	public void setRAS1(float val) {
		ras1 = val;
	}

	/** The relative allele strength for the second block */
	private float ras2;

	/**
	 * Gets the RAS2 value
	 * 
	 * @return The RAS2 value.
	 */
	public float getRAS2() {
		return ras2;
	}

	/**
	 * Sets the RAS2 value
	 * 
	 * @param val
	 *          The RAS2 value.
	 */
	public void setRAS2(float val) {
		ras2 = val;
	}

	/** The p-value associated with an AA call */
	private float pvalue_AA;

	/**
	 * Gets the AA p-value
	 * 
	 * @return The AA p-value.
	 */
	public float getPValue_AA() {
		return pvalue_AA;
	}

	/**
	 * Sets the AA p-value
	 * 
	 * @param val
	 *          The AA p-value.
	 */
	public void setPValue_AA(float val) {
		pvalue_AA = val;
	}

	/** The p-value associated with an AB call */
	private float pvalue_AB;

	/**
	 * Gets the AB p-value
	 * 
	 * @return The AB p-value.
	 */
	public float getPValue_AB() {
		return pvalue_AB;
	}

	/**
	 * Sets the AB p-value
	 * 
	 * @param val
	 *          The AB p-value.
	 */
	public void setPValue_AB(float val) {
		pvalue_AB = val;
	}

	/** The p-value associated with an BB call */
	private float pvalue_BB;

	/**
	 * Gets the BB p-value
	 * 
	 * @return The BB p-value.
	 */
	public float getPValue_BB() {
		return pvalue_BB;
	}

	/**
	 * Sets the BB p-value
	 * 
	 * @param val
	 *          The BB p-value.
	 */
	public void setPValue_BB(float val) {
		pvalue_BB = val;
	}

	/** The p-value associated with an no call call */
	private float pvalue_NoCall;

	/**
	 * Gets the NoCall p-value
	 * 
	 * @return The NoCall p-value.
	 */
	public float getPValue_NoCall() {
		return pvalue_NoCall;
	}

	/**
	 * Sets the NoCall p-value
	 * 
	 * @param val
	 *          The NoCall p-value.
	 */
	public void setPValue_NoCall(float val) {
		pvalue_NoCall = val;
	}

	/**
	 * Returns a string representation of the allele call.
	 * 
	 * @return The allele call as a string.
	 */
	public String getAlleleCallString() {
		switch (alleleCall) {
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
	public GenotypeProbeSetResults() {
		alleleCall = ALLELE_NO_CALL;
		confidence = 0.0f;
		ras1 = 0.0f;
		ras2 = 0.0f;
		pvalue_AA = 0.0f;
		pvalue_AB = 0.0f;
		pvalue_BB = 0.0f;
		pvalue_NoCall = 0.0f;

	}

}
