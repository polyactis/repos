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

/**
 * A type to hold the force call, the call that the algorithm would have made if the thresholds were relaxed.
 */
public class FusionForceCallType {

	/** The force call was made due to no signal threshold. */
	public static final byte NO_SIGNAL_THR_FORCE_CALL = 'N';

	/** The force call was made due to weak signal threshold. */
	public static final byte WEAK_SIGNAL_THR_FORCE_CALL = 'W';

	/** The force call was made due to saturation level. */
	public static final byte SATURATION_LEVEL_FORCE_CALL = 'S';

	/** The force call was made due to quality score threshold. */
	public static final byte QUALITY_SCORE_THR_FORCE_CALL = 'Q';

	/** The force call was made due to failed both trace and sequence profiles. */
	public static final byte TRACE_AND_SEQUENCE_PROFILES_FORCE_CALL = 'F';

	/** The force call was made due to base reliability threshold. */
	public static final byte RELIABILITY_THR_FORCE_CALL = 'B';

	/** The position (index) of the call. */
	private int position;

	/**
	 * Gets the position of the call.
	 * 
	 * @return The position.
	 */
	public int getPosition() {
		return position;
	}

	/**
	 * Sets the position of the call.
	 * 
	 * @param p
	 *          The position.
	 */
	public void setPosition(int p) {
		position = p;
	}

	/** The force call. */
	private byte call;

	/**
	 * Gets the call.
	 * 
	 * @return The call.
	 */
	public byte getCall() {
		return call;
	}

	/**
	 * Sets the call.
	 * 
	 * @param c
	 *          The call.
	 */
	public void setCall(byte c) {
		call = c;
	}

	/** The reason for the call. */
	private byte reason;

	/**
	 * Gets the reason of the call.
	 * 
	 * @return The reason.
	 */
	public byte getReason() {
		return reason;
	}

	/**
	 * Sets the reason of the call.
	 * 
	 * @param r
	 *          The reason.
	 */
	public void setReason(byte r) {
		reason = r;
	}

	/** Creates a new instance of FusionForceCallType */
	public FusionForceCallType() {
		position = 0;
		call = ' ';
		reason = ' ';
	}

	/** Creates a new instance of FusionForceCallType */
	public FusionForceCallType(FusionForceCallType f) {
		position = f.getPosition();
		call = f.getCall();
		reason = f.getReason();
	}

}
