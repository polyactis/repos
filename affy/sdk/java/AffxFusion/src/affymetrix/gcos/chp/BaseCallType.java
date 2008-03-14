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

/** Holds a call at a given position. */
public class BaseCallType {

	/** The position (index) of the call. */
	private int position;

	/**
	 * Gets the position (index) of the call.
	 * 
	 * @return The position
	 */
	public int getPosition() {
		return position;
	}

	/**
	 * Sets the position (index) of the call.
	 * 
	 * @param p
	 *          The position
	 */
	public void setPosition(int p) {
		position = p;
	}

	/** The call. */
	public char call;

	/**
	 * Gets the call.
	 * 
	 * @return The call
	 */
	public char getCall() {
		return call;
	}

	/**
	 * Sets the call.
	 * 
	 * @param c
	 *          The call
	 */
	public void setCall(char c) {
		call = c;
	}

	/** Creates a new instance of BaseCallType */
	public BaseCallType() {
		position = 0;
		call = ' ';
	}

	/** Creates a new instance of BaseCallType */
	public BaseCallType(BaseCallType b) {
		position = b.getPosition();
		call = b.getCall();
	}

}
