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

import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.data.CHPData;
import affymetrix.calvin.data.CHPReseqEntry;
import affymetrix.calvin.data.CHPReseqForceCall;
import affymetrix.calvin.data.CHPReseqOrigCall;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.gcos.chp.BaseCallType;
import affymetrix.gcos.chp.ForceCallType;
import affymetrix.gcos.chp.ResequencingResults;

/** Stores results for a resequencing analysis. */
public class FusionResequencingResults {

	/** The GCOS probe set results object. */
	private ResequencingResults gcosResult;

	/**
	 * Sets the GCOS probe set results object.
	 * 
	 * @param r
	 *          The GCOS probe set results object.
	 */
	public void setGCOSObject(ResequencingResults r) {
		gcosResult = r;
	}

	/** The calvin probe set results object. */
	private CHPData calvinResult;

	/**
	 * Sets the calvin object.
	 * 
	 * @param r
	 *          The calvin data file object.
	 */
	public void setCalvinObject(CHPData r) {
		calvinResult = r;
	}

	/** Clears the members. */
	public void clear() {
		gcosResult = null;
		calvinResult = null;
	}

	/** Creates a new instance of FusionResequencingResults */
	public FusionResequencingResults() {
		clear();
	}

	/**
	 * Gets the called bases.
	 * 
	 * @return The array of called bases.
	 */
	public List<Byte> getCalledBases() throws UnsignedOutOfLimitsException {
		if (gcosResult != null) {
			List<Byte> calledBases = new ArrayList<Byte>();
			int n = gcosResult.getCalledBasesSize();
			for (int i = 0; i < n; i++) {
				calledBases.add(i, (byte)gcosResult.getCalledBase(i));
			}
			return calledBases;
		}
		else if (calvinResult != null) {
			int n = calvinResult.getEntryCount();
			List<Byte> calledBases = new ArrayList<Byte>();
			CHPReseqEntry entry = new CHPReseqEntry();
			for (int i = 0; i < n; i++) {
				calvinResult.getEntry(i, entry);
				calledBases.add(i, entry.getCall());
			}
			return calledBases;
		}
		return null;
	}

	/**
	 * Gets the called base at the given index.
	 * 
	 * @param index
	 *          The index to the called bases array.
	 * @return The called base.
	 */
	public byte getCalledBase(int index) throws UnsignedOutOfLimitsException {
		if (gcosResult != null) {
			return (byte)gcosResult.getCalledBase(index);
		}
		else if (calvinResult != null) {
			CHPReseqEntry entry = new CHPReseqEntry();
			calvinResult.getEntry(index, entry);
			return entry.getCall();
		}
		return ' ';
	}

	/**
	 * Gets the size of the called bases array.
	 * 
	 * @return The size of the called bases array.
	 */
	public int getCalledBasesSize() {
		if (gcosResult != null) {
			return gcosResult.getCalledBasesSize();
		}
		else if (calvinResult != null) {
			return calvinResult.getEntryCount();
		}
		return 0;
	}

	/**
	 * Gets the scores.
	 * 
	 * @return The array of scores.
	 */
	public List<Float> getScores() throws UnsignedOutOfLimitsException {
		if (gcosResult != null) {
			List<Float> scores = new ArrayList<Float>();
			int n = gcosResult.getScoresSize();
			for (int i = 0; i < n; i++) {
				scores.add(i, gcosResult.getScore(i));
			}
			return scores;
		}
		else if (calvinResult != null) {
			int n = calvinResult.getEntryCount();
			List<Float> scores = new ArrayList<Float>();
			CHPReseqEntry entry = new CHPReseqEntry();
			for (int i = 0; i < n; i++) {
				calvinResult.getEntry(i, entry);
				scores.add(i, entry.getScore());
			}
			return scores;
		}
		return null;
	}

	/**
	 * Gets the score at the given index.
	 * 
	 * @param index
	 *          The index to the scores array.
	 * @return The score.
	 */
	public float getScore(int index) throws UnsignedOutOfLimitsException {
		if (gcosResult != null) {
			return gcosResult.getScore(index);
		}
		else if (calvinResult != null) {
			CHPReseqEntry entry = new CHPReseqEntry();
			calvinResult.getEntry(index, entry);
			return entry.getScore();
		}
		return 0.0f;
	}

	/**
	 * Gets the size of the scores array.
	 * 
	 * @return The size of the scores array.
	 */
	public int getScoresSize() {
		if (gcosResult != null) {
			return gcosResult.getScoresSize();
		}
		else if (calvinResult != null) {
			return calvinResult.getEntryCount();
		}
		return 0;
	}

	/**
	 * Gets the force calls.
	 * 
	 * @return The array of force calls.
	 */
	public List<FusionForceCallType> getForceCalls() throws UnsignedOutOfLimitsException {
		if (gcosResult != null) {
			List<FusionForceCallType> force = new ArrayList<FusionForceCallType>();
			int n = gcosResult.getForceCallsSize();
			for (int i = 0; i < n; i++) {
				ForceCallType f = gcosResult.getForceCall(i);
				FusionForceCallType ff = new FusionForceCallType();
				ff.setCall(f.getCall());
				ff.setPosition(f.getPosition());
				ff.setReason(f.getReason());
				force.set(i, ff);
			}
			return force;
		}
		else if (calvinResult != null) {
			List<FusionForceCallType> force = new ArrayList<FusionForceCallType>();
			int n = calvinResult.getForceCnt();
			CHPReseqForceCall f = new CHPReseqForceCall();
			for (int i = 0; i < n; i++) {
				calvinResult.getForceCall(i, f);
				FusionForceCallType ff = new FusionForceCallType();
				ff.setCall(f.getCall());
				ff.setPosition(f.getPosition());
				ff.setReason(f.getReason());
				force.set(i, ff);
			}
			return force;
		}
		return null;
	}

	/**
	 * Gets the force call at the given index.
	 * 
	 * @param index
	 *          The index to the force calls array.
	 * @return The force call.
	 */
	public FusionForceCallType getForceCall(int index) throws UnsignedOutOfLimitsException {
		if (gcosResult != null) {
			ForceCallType gcosForce = gcosResult.getForceCall(index);
			FusionForceCallType force = new FusionForceCallType();
			force.setCall(gcosForce.getCall());
			force.setPosition(gcosForce.getPosition());
			force.setReason(gcosForce.getReason());
			return force;
		}
		else if (calvinResult != null) {
			CHPReseqForceCall f = new CHPReseqForceCall();
			calvinResult.getForceCall(index, f);
			FusionForceCallType force = new FusionForceCallType();
			force.setCall(f.getCall());
			force.setPosition(f.getPosition());
			force.setReason(f.getReason());
			return force;
		}
		return null;
	}

	/**
	 * Gets the size of the force calls array.
	 * 
	 * @return The size of the force calls array.
	 */
	public int getForceCallsSize() {
		if (gcosResult != null) {
			return gcosResult.getForceCallsSize();
		}
		else if (calvinResult != null) {
			return calvinResult.getForceCnt();
		}
		return 0;
	}

	/**
	 * Gets the original called bases.
	 * 
	 * @return The array of original calls.
	 */
	public List<FusionBaseCallType> getOrigCalls() throws UnsignedOutOfLimitsException {
		if (gcosResult != null) {
			List<FusionBaseCallType> bases = new ArrayList<FusionBaseCallType>();
			int n = gcosResult.getOrigCallsSize();
			for (int i = 0; i < n; i++) {
				BaseCallType gcosBase = gcosResult.getOrigCall(i);
				FusionBaseCallType base = new FusionBaseCallType();
				base.setCall((byte)gcosBase.getCall());
				base.setPosition(gcosBase.getPosition());
				bases.set(i, base);
			}
			return bases;
		}
		else if (calvinResult != null) {
			List<FusionBaseCallType> bases = new ArrayList<FusionBaseCallType>();
			int n = calvinResult.getOrigCnt();
			CHPReseqOrigCall calvinBase = new CHPReseqOrigCall();
			for (int i = 0; i < n; i++) {
				calvinResult.getOrigCall(i, calvinBase);
				FusionBaseCallType base = new FusionBaseCallType();
				base.setCall(calvinBase.getCall());
				base.setPosition(calvinBase.getPosition());
				bases.set(i, base);
			}
			return bases;
		}
		return null;
	}

	/**
	 * Gets the original called base at the given index.
	 * 
	 * @param index
	 *          The index to the original calls array.
	 * @return The original call.
	 */
	public FusionBaseCallType getOrigCall(int index) throws UnsignedOutOfLimitsException {
		if (gcosResult != null) {
			BaseCallType gcosBase = gcosResult.getOrigCall(index);
			FusionBaseCallType base = new FusionBaseCallType();
			base.setCall((byte)gcosBase.getCall());
			base.setPosition(gcosBase.getPosition());
			return base;
		}
		else if (calvinResult != null) {
			CHPReseqOrigCall calvinBase = new CHPReseqOrigCall();
			calvinResult.getOrigCall(index, calvinBase);
			FusionBaseCallType base = new FusionBaseCallType();
			base.setCall(calvinBase.getCall());
			base.setPosition(calvinBase.getPosition());
			return base;
		}
		return null;
	}

	/**
	 * Gets the size of the original calls array.
	 * 
	 * @return The size of the original calls array.
	 */
	public int getOrigCallsSize() {
		if (gcosResult != null) {
			return gcosResult.getOrigCallsSize();
		}
		else if (calvinResult != null) {
			return calvinResult.getOrigCnt();
		}
		return 0;
	}
}
