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

package affymetrix.fusion.chp;

import affymetrix.calvin.data.CHPExpressionEntry;
import affymetrix.gcos.chp.ExpressionProbeSetResults;
import affymetrix.portability.UByte;
import affymetrix.portability.UShort;

/** Expression analysis probe set results for the MAS5 algorithm */
public class FusionExpressionProbeSetResults {

	/** Present call for expression analysis */
	public static final byte ABS_PRESENT_CALL = 0;

	/** Marginal call for expression analysis */
	public static final byte ABS_MARGINAL_CALL = 1;

	/** Absent call for expression analysis */
	public static final byte ABS_ABSENT_CALL = 2;

	/** No call call for expression analysis */
	public static final byte ABS_NO_CALL = 3;

	/** Increase call for expression comparison analysis */
	public static final byte COMP_INCREASE_CALL = 1;

	/** Decrease call for expression comparison analysis */
	public static final byte COMP_DECREASE_CALL = 2;

	/** Moderate increase call for expression comparison analysis */
	public static final byte COMP_MOD_INCREASE_CALL = 3;

	/** Moderate decrease call for expression comparison analysis */
	public static final byte COMP_MOD_DECREASE_CALL = 4;

	/** No change call for expression comparison analysis */
	public static final byte COMP_NO_CHANGE_CALL = 5;

	/** No call call for expression comparison analysis */
	public static final byte COMP_NO_CALL = 6;

	/** Creates a new instance of FusionExpressionProbeSetResults */
	public FusionExpressionProbeSetResults() {
		clear();
	}

	/** The GCOS probe set results object. */
	private ExpressionProbeSetResults gcosResult;

	/**
	 * Sets the GCOS probe set results object.
	 * 
	 * @param r
	 *          The GCOS probe set results object.
	 */
	public void setGCOSObject(ExpressionProbeSetResults r) {
		gcosResult = r;
	}

	/** The calvin probe set results object. */
	private CHPExpressionEntry calvinResult;

	/**
	 * Sets the calvin object.
	 * 
	 * @param r
	 *          The calvin probe set object.
	 */
	public void setCalvinObject(CHPExpressionEntry r) {
		calvinResult = r;
	}

	/** Clears the members. */
	public void clear() {
		gcosResult = null;
		calvinResult = null;
	}

	/**
	 * Gets detection pvalue return Detection pvalue.
	 */
	public float getDetectionPValue() {
		if (gcosResult != null) {
			return gcosResult.getDetectionPValue();
		}
		else if (calvinResult != null) {
			return calvinResult.getDetectionPValue();
		}
		return 0.0f;
	}

	/**
	 * Gets the signal. return Signal.
	 */
	public float getSignal() {
		if (gcosResult != null) {
			return gcosResult.getSignal();
		}
		else if (calvinResult != null) {
			return calvinResult.getSignal();
		}
		return 0.0f;
	}

	/**
	 * gets the number of pairs. return Number of pairs.
	 */
	public UShort getNumPairs() {
		if (gcosResult != null) {
			return gcosResult.getNumPairs();
		}
		else if (calvinResult != null) {
			return calvinResult.getNumPairs();
		}
		return UShort.ZERO;
	}

	/**
	 * Gets the number used pairs. return Number of used pairs.
	 */
	public UShort getNumUsedPairs() {
		if (gcosResult != null) {
			return gcosResult.getNumUsedPairs();
		}
		else if (calvinResult != null) {
			return calvinResult.getNumPairsUsed();
		}
		return UShort.ZERO;
	}

	/**
	 * Gets the detection. return Detection.
	 */
	public UByte getDetection() {
		if (gcosResult != null) {
			return gcosResult.getDetection();
		}
		else if (calvinResult != null) {
			return calvinResult.getDetection();
		}
		return UByte.ZERO;
	}

	/**
	 * Determines if comp results exsist. return True if comp results exsist.
	 */
	public boolean hasCompResults() {
		if (gcosResult != null) {
			return gcosResult.getHasCompResults();
		}
		else if (calvinResult != null) {
			return calvinResult.getHasComparisonData();
		}
		return false;
	}

	/**
	 * Gets the change pvalue. return Cahnge pvalue.
	 */
	public float getChangePValue() {
		if (gcosResult != null) {
			return gcosResult.getChangePValue();
		}
		else if (calvinResult != null) {
			return calvinResult.getChangePValue();
		}
		return 0.0f;
	}

	/**
	 * Gets the signal log ratio. return Signal log ratio.
	 */
	public float getSignalLogRatio() {
		if (gcosResult != null) {
			return gcosResult.getSignalLogRatio();
		}
		else if (calvinResult != null) {
			return calvinResult.getSigLogRatio();
		}
		return 0.0f;
	}

	/**
	 * Gets the signal log ratio low return Signal log ratio low.
	 */
	public float getSignalLogRatioLow() {
		if (gcosResult != null) {
			return gcosResult.getSignalLogRatioLow();
		}
		else if (calvinResult != null) {
			return calvinResult.getSigLogRatioLo();
		}
		return 0.0f;
	}

	/**
	 * Gets the signal log ratio high return Signal log ratio high.
	 */
	public float getSignalLogRatioHigh() {
		if (gcosResult != null) {
			return gcosResult.getSignalLogRatioHigh();
		}
		else if (calvinResult != null) {
			return calvinResult.getSigLogRatioHi();
		}
		return 0.0f;
	}

	/**
	 * Gets the number of common pairs return Number of common pairs.
	 */
	public UShort getNumCommonPairs() {
		if (gcosResult != null) {
			return gcosResult.getNumCommonPairs();
		}
		else if (calvinResult != null) {
			return calvinResult.getCommonPairs();
		}
		return new UShort();
	}

	/**
	 * Gets the change. return Change.
	 */
	public UByte getChange() {
		if (gcosResult != null) {
			return gcosResult.getChange();
		}
		else if (calvinResult != null) {
			return calvinResult.getChange();
		}
		return new UByte();
	}

	/**
	 * Returns a string representation of the detection call.
	 * 
	 * @return The detection call
	 */
	public String getDetectionString() {
		switch (getDetection().toShort()) {
		case ABS_PRESENT_CALL:
			return "P";

		case ABS_MARGINAL_CALL:
			return "M";

		case ABS_ABSENT_CALL:
			return "A";

		case ABS_NO_CALL:
			return "No Call";
		}
		return null;
	}

	/**
	 * Returns a string representation of the change call.
	 * 
	 * @return The change call
	 */
	public String getChangeString() {
		switch (getChange().toShort()) {
		case COMP_INCREASE_CALL:
			return "I";

		case COMP_DECREASE_CALL:
			return "D";

		case COMP_MOD_INCREASE_CALL:
			return "MI";

		case COMP_MOD_DECREASE_CALL:
			return "MD";

		case COMP_NO_CHANGE_CALL:
			return "NC";

		case COMP_NO_CALL:
			return "No Call";
		}
		return null;
	}
}
