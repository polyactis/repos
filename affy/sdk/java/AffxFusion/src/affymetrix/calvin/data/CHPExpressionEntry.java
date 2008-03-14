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

import affymetrix.portability.UByte;
import affymetrix.portability.UShort;

/** This class stores the expression probe set analysis results. */
public class CHPExpressionEntry {

	/** The probe set name. */
	private String probeSetName;

	/** The detection call. */
	private UByte detection;

	/** The detection p-value. */
	private float detectionPValue;

	/** The signal. */
	private float signal;

	/** The number of probe pairs in the probe set. */
	private UShort numPairs;

	/** The number of probe pairs used in the analysis. */
	private UShort numPairsUsed;

	/** A flag to indicate if comparison data exists. */
	private boolean hasComparisonData;

	/** The change call. */
	private UByte change;

	/** The change p-value. */
	private float changePValue;

	/** The signal log ratio. */
	private float sigLogRatio;

	/** The signal log ratio low. */
	private float sigLogRatioLo;

	/** The signal log ratio high. */
	private float sigLogRatioHi;

	/** The number of probe pairs in common between control and experiment. */
	private UShort commonPairs;

	/** Constructor. */
	public CHPExpressionEntry() {
		detection = UByte.ZERO;
		detectionPValue = 0;
		signal = 0;
		numPairs = UShort.ZERO;
		numPairsUsed = UShort.ZERO;
		hasComparisonData = false;
		change = UByte.ZERO;
		changePValue = 0;
		sigLogRatio = 0;
		sigLogRatioLo = 0;
		sigLogRatioHi = 0;
		commonPairs = UShort.ZERO;
	}

	/**
	 * Constructor with absolute and comparison data.
	 * 
	 * @param psName
	 *          The probe set name
	 * @param detect
	 *          The detection
	 * @param detectPValue
	 *          The detection p-value
	 * @param sig
	 *          The signal
	 * @param nPairs
	 *          The number of probe pairs in the set.
	 * @param nPairsUsed
	 *          The number of probe pairs used in the analysis
	 * @param compData
	 *          Flag indicating if comp data exists.
	 * @param chg
	 *          The change call
	 * @param chgPValue
	 *          The change p-value
	 * @param sLogRatio
	 *          The signal log ratio
	 * @param sLogRatioLo
	 *          The signal log ratio low
	 * @param sLogRatioHi
	 *          The signal log ratio high
	 * @param commonPrs
	 *          The number of probe pairs in common between control and experiment
	 */
	public CHPExpressionEntry(String psName, UByte detect, float detectPValue, float sig, UShort nPairs,
			UShort nPairsUsed, boolean compData, UByte chg, float chgPValue, float sLogRatio, float sLogRatioLo,
			float sLogRatioHi, UShort commonPrs) {
		setProbeSetName(psName);
		setDetection(detect);
		setDetectionPValue(detectPValue);
		setSignal(sig);
		setNumPairs(nPairs);
		setNumPairsUsed(nPairsUsed);
		setHasComparisonData(compData);
		setChange(chg);
		setChangePValue(chgPValue);
		setSigLogRatio(sLogRatio);
		setSigLogRatioLo(sLogRatioLo);
		setSigLogRatioHi(sLogRatioHi);
		setCommonPairs(commonPrs);
	}

	/**
	 * Constructor with absolute data only.
	 * 
	 * @param psName
	 *          The probe set name
	 * @param detect
	 *          The detection
	 * @param detectPValue
	 *          The detection p-value
	 * @param sig
	 *          The signal
	 * @param nPairs
	 *          The number of probe pairs in the set.
	 * @param nPairsUsed
	 *          The number of probe pairs used in the analysis
	 */
	public CHPExpressionEntry(String psName, UByte detect, float detectPValue, float sig, UShort nPairs, UShort nPairsUsed) {
		this(psName, detect, detectPValue, sig, nPairs, nPairsUsed, false, new UByte(), 0f, 0f, 0f, 0f, new UShort());
	}

	/**
	 * Constructor with expression probe set result.
	 * 
	 * @param e
	 *          The expression result.
	 */
	public CHPExpressionEntry(CHPExpressionEntry e) {
		this(e.getProbeSetName(), e.getDetection(), e.getDetectionPValue(), e.getSignal(), e.getNumPairs(), e
				.getNumPairsUsed(), e.getHasComparisonData(), e.getChange(), e.getChangePValue(), e.getSigLogRatio(), e
				.getSigLogRatioLo(), e.getSigLogRatioHi(), e.getCommonPairs());
	}

	/** Clears the members. */
	public void clear() {
		probeSetName = null;
		detection = UByte.ZERO;
		detectionPValue = 0;
		signal = 0;
		numPairs = UShort.ZERO;
		numPairsUsed = UShort.ZERO;
		hasComparisonData = false;
		change = UByte.ZERO;
		changePValue = 0;
		sigLogRatio = 0;
		sigLogRatioLo = 0;
		sigLogRatioHi = 0;
		commonPairs = UShort.ZERO;
	}

	/** Gets the probe set name. */
	public String getProbeSetName() {
		return probeSetName;
	}

	/** Get the detection call. */
	public UByte getDetection() {
		return detection;
	}

	/** Gets the detection p-value. */
	public float getDetectionPValue() {
		return detectionPValue;
	}

	/** Gets the signal value. */
	public float getSignal() {
		return signal;
	}

	/** Gets the number of probe pairs in the set. */
	public UShort getNumPairs() {
		return numPairs;
	}

	/** Gets the number of probe pairs used in the analysis. */
	public UShort getNumPairsUsed() {
		return numPairsUsed;
	}

	/** A flag to indicate if comparison data exists. */
	public boolean getHasComparisonData() {
		return hasComparisonData;
	}

	/** Gets the change call. */
	public UByte getChange() {
		return change;
	}

	/** Gets the change p-value. */
	public float getChangePValue() {
		return changePValue;
	}

	/** Gets the signal log ratio. */
	public float getSigLogRatio() {
		return sigLogRatio;
	}

	/** Gets the signal log ratio low. */
	public float getSigLogRatioLo() {
		return sigLogRatioLo;
	}

	/** Gets the signal log ratio high. */
	public float getSigLogRatioHi() {
		return sigLogRatioHi;
	}

	/** Gets the number of probe pairs in common between control and experiment. */
	public UShort getCommonPairs() {
		return commonPairs;
	}

	/** Sets the probe set name. */
	public void setProbeSetName(String p) {
		probeSetName = p;
	}

	/** Sets the detection value. */
	public void setDetection(UByte p) {
		detection = p;
	}

	/** Sets the detection p-value. */
	public void setDetectionPValue(float p) {
		detectionPValue = p;
	}

	/** Sets the signal value. */
	public void setSignal(float p) {
		signal = p;
	}

	/** Sets the number of pairs in the probe sets. */
	public void setNumPairs(UShort p) {
		numPairs = p;
	}

	/** Sets the number of pairs used in the analysis. */
	public void setNumPairsUsed(UShort p) {
		numPairsUsed = p;
	}

	/** Sets the flag to indicate if comparison data exists. */
	public void setHasComparisonData(boolean b) {
		hasComparisonData = b;
	}

	/** Sets the change call. */
	public void setChange(UByte p) {
		change = p;
	}

	/** Sets the change p-value. */
	public void setChangePValue(float p) {
		changePValue = p;
	}

	/** Sets the signal log ratio. */
	public void setSigLogRatio(float p) {
		sigLogRatio = p;
	}

	/** Sets the signal log ratio low. */
	public void setSigLogRatioLo(float p) {
		sigLogRatioLo = p;
	}

	/** Sets the signal log ratio high. */
	public void setSigLogRatioHi(float p) {
		sigLogRatioHi = p;
	}

	/** Sets the number of probe pairs in common between control and experiment. */
	public void setCommonPairs(UShort p) {
		commonPairs = p;
	}

}
