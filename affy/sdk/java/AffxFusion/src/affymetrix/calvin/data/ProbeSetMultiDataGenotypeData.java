////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
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
////////////////////////////////////////////////////////////////

package affymetrix.calvin.data;

/** Stores data for a genotype result of a probe set. */
public class ProbeSetMultiDataGenotypeData extends ProbeSetMultiDataBase {

	/* ! The AA call */
	public static final byte SNP_AA_CALL = 6;

	/* ! The BB call */
	public static final byte SNP_BB_CALL = 7;

	/* ! The AB call */
	public static final byte SNP_AB_CALL = 8;

	/* ! The no call allele call */
	public static final byte SNP_NO_CALL = 11;

	/** The name of the probe set. */
	private String name;

	/** The call. */
	private byte call;

	/** The confidence of the call. */
	private float confidence;

	/** Creates a new instance of ProbeSetMultiDataGenotypeData */
	public ProbeSetMultiDataGenotypeData() {
		name = null;
		call = 0;
		confidence = 0.0f;
	}

	/** The name of the probe set. */
	public String getName() {
		return name;
	}

	/** Set the name */
	public void setName(String n) {
		name = n;
	}

	/** Get the call */
	public byte getCall() {
		return call;
	}

	/** Set the call. */
	public void setCall(byte c) {
		call = c;
	}

	/** Get the confidence. */
	public float getConfidence() {
		return confidence;
	}

	/** Set the confidence */
	public void setConfidence(float c) {
		confidence = c;
	}

	/**
	 * Convert a string representation of a genotype call to a numeric representation.
	 * 
	 * @param call
	 *          The call value.
	 * @return A numeric representation of the call value.
	 */
	public static byte genotypeCallFromString(String call) {
		if ((call.equals("A") == true) || (call.equals("AA") == true)) {
			return SNP_AA_CALL;
		}
		else if ((call.equals("B") == true) || (call.equals("BB") == true)) {
			return SNP_BB_CALL;
		}
		else if (call.equals("AB") == true) {
			return SNP_AB_CALL;
		}
		else {
			return SNP_NO_CALL;
		}
	}

	/**
	 * Convert a numeric representation of a genotype call to a string representation.
	 * 
	 * @param call
	 *          The call value.
	 * @return A string representation of the call value.
	 */
	public static String genotypeCallToString(byte call) {
		switch (call) {
		case SNP_AA_CALL:
			return "A";

		case SNP_AB_CALL:
			return "AB";

		case SNP_BB_CALL:
			return "BB";

		default:
			return "No Call";
		}
	}

}
