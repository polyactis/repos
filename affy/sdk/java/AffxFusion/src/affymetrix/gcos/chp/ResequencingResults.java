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

package affymetrix.gcos.chp;

import java.util.ArrayList;
import java.util.List;

/** Stores results for a resequencing analysis. */
public class ResequencingResults {

	/** The called bases. */
	private List<Character> calledBases = new ArrayList<Character>();

	/** Base call scores. */
	private List<Float> scores = new ArrayList<Float>();

	/** An array of force calls - base calls the algorithm would have made if the thresholds were removed. */
	private List<ForceCallType> forceCalls = new ArrayList<ForceCallType>();

	/**
	 * An array of original calls. The calledBases contained the results of the algorithm and user edits. If a user edits
	 * a base the original algorithm called base is stored in this vector.
	 */
	private List<BaseCallType> origCalls = new ArrayList<BaseCallType>();

	/** Creates a new instance of ResequencingResults */
	public ResequencingResults() {
	}

	/** Clears the members. */
	public void clear() {
		calledBases.clear();
		scores.clear();
		forceCalls.clear();
		origCalls.clear();
	}

	/**
	 * Gets the called bases.
	 * 
	 * @return The array of called bases.
	 */
	public List<Character> getCalledBases() {
		return calledBases;
	}

	/**
	 * Gets the called base at the given index.
	 * 
	 * @param index
	 *          The index to the called bases array.
	 * @return The called base.
	 */
	public char getCalledBase(int index) {
		return calledBases.get(index);
	}

	/**
	 * Gets the size of the called bases array.
	 * 
	 * @return The size of the called bases array.
	 */
	public int getCalledBasesSize() {
		return calledBases.size();
	}

	public void clearCalledBases() {
		calledBases.clear();
	}

	/**
	 * Sets the called base.
	 * 
	 * @param index
	 *          The index to the array.
	 * @param call
	 *          The call.
	 */
	public void addCalledBase(char call) {
		calledBases.add(call);
	}

	/**
	 * Gets the scores.
	 * 
	 * @return The array of scores.
	 */
	public List<Float> getScores() {
		return scores;
	}

	/**
	 * Gets the score at the given index.
	 * 
	 * @param index
	 *          The index to the scores array.
	 * @return The score.
	 */
	public float getScore(int index) {
		return scores.get(index);
	}

	/**
	 * Gets the size of the scores array.
	 * 
	 * @return The size of the scores array.
	 */
	public int getScoresSize() {
		return scores.size();
	}

	public void clearScores() {
		scores.clear();
	}

	/**
	 * Sets the score.
	 * 
	 * @param index
	 *          The index to the array.
	 * @param score
	 *          The score.
	 */
	public void addScore(float score) {
		scores.add(score);
	}

	/**
	 * Gets the force calls.
	 * 
	 * @return The array of force calls.
	 */
	public List<ForceCallType> getForceCalls() {
		return forceCalls;
	}

	/**
	 * Gets the force call at the given index.
	 * 
	 * @param index
	 *          The index to the force calls array.
	 * @return The force call.
	 */
	public ForceCallType getForceCall(int index) {
		return forceCalls.get(index);
	}

	/**
	 * Gets the size of the force calls array.
	 * 
	 * @return The size of the force calls array.
	 */
	public int getForceCallsSize() {
		return forceCalls.size();
	}

	/**
	 * Resizes the force calls array.
	 * 
	 * @param size
	 *          The size of the array.
	 */
	public void clearForceCalls() {
		forceCalls.clear();
	}

	/**
	 * Sets the force call.
	 * 
	 * @param index
	 *          The index to the array.
	 * @param call
	 *          The force call.
	 */
	public void addForceCall(ForceCallType call) {
		forceCalls.add(call);
	}

	/**
	 * Gets the original called bases.
	 * 
	 * @return The array of original calls.
	 */
	public List<BaseCallType> getOrigCalls() {
		return origCalls;
	}

	/**
	 * Gets the original called base at the given index.
	 * 
	 * @param index
	 *          The index to the original calls array.
	 * @return The original call.
	 */
	public BaseCallType getOrigCall(int index) {
		return origCalls.get(index);
	}

	/**
	 * Gets the size of the original calls array.
	 * 
	 * @return The size of the original calls array.
	 */
	public int getOrigCallsSize() {
		return origCalls.size();
	}

	/**
	 * Resizes the original calls array.
	 * 
	 * @param size
	 *          The size of the array.
	 */
	public void clearOrigCalls() {
		origCalls.clear();
	}

	/**
	 * Sets the original call.
	 * 
	 * @param index
	 *          The index to the array.
	 * @param call
	 *          The original call.
	 */
	public void addOrigCall(BaseCallType call) {
		origCalls.add(call);
	}

}
