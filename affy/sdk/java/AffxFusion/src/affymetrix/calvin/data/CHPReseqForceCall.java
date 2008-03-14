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

/**
 * A structure to hold a force call, its position and reason.
 * 
 * A force call is the call the algorithm would have made if the thresholds were not applied.
 */
public class CHPReseqForceCall {

	public CHPReseqForceCall() {
		position = 0;
		call = ' ';
		reason = ' ';
	}

	/** The position (index) of the call. */
	private int position;

	public int getPosition() {
		return position;
	}

	public void setPosition(int p) {
		position = p;
	}

	/** The call at the given position. */
	private byte call;

	public byte getCall() {
		return call;
	}

	public void setCall(byte c) {
		call = c;
	}

	/** The reason for the call. */
	private byte reason;

	public byte getReason() {
		return reason;
	}

	public void setReason(byte r) {
		reason = r;
	}

}
