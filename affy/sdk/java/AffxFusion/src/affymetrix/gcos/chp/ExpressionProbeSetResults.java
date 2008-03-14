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

import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.portability.UByte;
import affymetrix.portability.UShort;

/** Defines the results for an expression probe set. */
public class ExpressionProbeSetResults {

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

	/** The detection p-value */
	private float detectionPValue;

	/**
	 * Gets the detection p-value
	 * 
	 * @return The detection p-value.
	 */
	public float getDetectionPValue() {
		return detectionPValue;
	}

	/**
	 * Sets the detection p-value
	 * 
	 * @param p
	 *          The detection p-value.
	 */
	public void setDetectionPValue(float p) {
		detectionPValue = p;
	}

	/** The signal value */
	private float signal;

	/**
	 * Gets the signal
	 * 
	 * @return The signal.
	 */
	public float getSignal() {
		return signal;
	}

	/**
	 * Sets the signal
	 * 
	 * @param s
	 *          The signal.
	 */
	public void setSignal(float s) {
		signal = s;
	}

	/** The number of probe pairs in the set */
	private UShort numPairs;

	/**
	 * Gets the number of probe pairs in the set.
	 * 
	 * @return The number of pairs.
	 */
	public UShort getNumPairs() {
		return numPairs;
	}

	/**
	 * Sets the number of probe pairs in the set.
	 * 
	 * @param n
	 *          The number of pairs.
	 */
	public void setNumPairs(UShort n) {
		numPairs = n;
	}

	/** The number of probe pairs used to calculate the signal value */
	private UShort numUsedPairs;

	/**
	 * Gets the number of probe pairs used to calculate the signal.
	 * 
	 * @return The number of pairs used to calculate the signal.
	 */
	public UShort getNumUsedPairs() {
		return numUsedPairs;
	}

	/**
	 * Sets the number of probe pairs used to calculate the signal.
	 * 
	 * @param n
	 *          The number of pairs used to calculate the signal.
	 */
	public void setNumUsedPairs(UShort n) {
		numUsedPairs = n;
	}

	/** The detection call */
	private UByte detection;

	/**
	 * Gets the detection.
	 * 
	 * @return The detection.
	 */
	public UByte getDetection() {
		return detection;
	}

	/**
	 * Sets the detection.
	 * 
	 * @param d
	 *          The detection.
	 */
	public void setDetection(UByte d) {
		detection = d;
	}

	/** Flag indicating that comparison results exist */
	private boolean hasCompResults;

	/**
	 * Gets a flag that indicates if comparison data exists.
	 * 
	 * @return The flag that indicates if comparison data exists.
	 */
	public boolean getHasCompResults() {
		return hasCompResults;
	}

	/**
	 * Sets a flag that indicates if comparison data exists.
	 * 
	 * @param b
	 *          The flag that indicates if comparison data exists.
	 */
	public void setHasCompResults(boolean b) {
		hasCompResults = b;
	}

	/** The change p-value */
	private float changePValue;

	/**
	 * Gets the change p-value
	 * 
	 * @return The change p-value
	 */
	public float getChangePValue() {
		return changePValue;
	}

	/**
	 * Set the change p-value
	 * 
	 * @param p
	 *          The change p-value.
	 */
	public void setChangePValue(float p) {
		changePValue = p;
	}

	/** The signal log ratio */
	private float signalLogRatio;

	/**
	 * Gets the signal log ratio.
	 * 
	 * @return The signal log ratio.
	 */
	public float getSignalLogRatio() {
		return signalLogRatio;
	}

	/**
	 * Set the signal log ratio.
	 * 
	 * @param s
	 *          The signal log ratio.
	 */
	public void setSignalLogRatio(float s) {
		signalLogRatio = s;
	}

	/** The signal log ratio low value */
	private float signalLogRatioLow;

	/**
	 * Gets the signal log ratio low.
	 * 
	 * @return The signal log ratio low.
	 */
	public float getSignalLogRatioLow() {
		return signalLogRatioLow;
	}

	/**
	 * Set the signal log ratio low.
	 * 
	 * @param s
	 *          The signal log ratio low.
	 */
	public void setSignalLogRatioLow(float s) {
		signalLogRatioLow = s;
	}

	/** The signal log ratio high value */
	private float signalLogRatioHigh;

	/**
	 * Gets the signal log ratio high.
	 * 
	 * @return The signal log ratio high.
	 */
	public float getSignalLogRatioHigh() {
		return signalLogRatioHigh;
	}

	/**
	 * Set the signal log ratio high.
	 * 
	 * @param s
	 *          The signal log ratio high.
	 */
	public void setSignalLogRatioHigh(float s) {
		signalLogRatioHigh = s;
	}

	/** The number of probe pairs in common between the experiment and baseline data */
	private UShort numCommonPairs;

	/**
	 * Gets the number of probe pairs in common between the experiment and baseline data.
	 * 
	 * @return The number of probe pairs in common between the experiment and baseline data.
	 */
	public UShort getNumCommonPairs() {
		return numCommonPairs;
	}

	/**
	 * Set the number of probe pairs in common between the experiment and baseline data.
	 * 
	 * @param n
	 *          The number of probe pairs in common between the experiment and baseline data.
	 */
	public void setNumCommonPairs(UShort n) {
		numCommonPairs = n;
	}

	/** The change call */
	private UByte change;

	/**
	 * Gets the change call.
	 * 
	 * @return The change call.
	 */
	public UByte getChange() {
		return change;
	}

	/**
	 * Set the change call.
	 * 
	 * @param c
	 *          The the change call.
	 */
	public void setChange(UByte c) {
		change = c;
	}

	/**
	 * Returns a string representation of the detection call.
	 * 
	 * @return The detection call
	 */
	public String getDetectionString() {
		switch (detection.toShort()) {
		case ABS_PRESENT_CALL:
			return "P";

		case ABS_MARGINAL_CALL:
			return "M";

		case ABS_ABSENT_CALL:
			return "A";

		case ABS_NO_CALL:
			return "No Call";

		default:
			return "";
		}
	}

	/**
	 * Returns a string representation of the change call.
	 * 
	 * @return The change call
	 */
	public String getChangeString() {
		switch (change.toShort()) {
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

		default:
			return "";
		}
	}

	/** Creates a new instance of ExpressionProbeSetResults */
	public ExpressionProbeSetResults() {
		try {
			detectionPValue = 0.0f;
			signal = 0.0f;
			numPairs = UShort.ZERO;
			numUsedPairs = UShort.ZERO;
			detection = new UByte(ABS_NO_CALL);
			hasCompResults = false;
			changePValue = 0.0f;
			signalLogRatio = 0.0f;
			signalLogRatioLow = 0.0f;
			signalLogRatioHigh = 0.0f;
			numCommonPairs = UShort.ZERO;
			change = new UByte(COMP_NO_CALL);
		}
		catch (UnsignedOutOfLimitsException e) {
			e.printStackTrace();
		}
	}
}
