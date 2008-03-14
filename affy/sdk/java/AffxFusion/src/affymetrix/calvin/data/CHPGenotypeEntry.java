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

/** This class stores a zone's background value */
public class CHPGenotypeEntry {

	private String probeSetName;

	private byte call;

	private float confidence;

	private float RAS1;

	private float RAS2;

	private float aaCall;

	private float abCall;

	private float bbCall;

	private float noCall;

	public CHPGenotypeEntry() {
		clear();
	}

	public CHPGenotypeEntry(String psName, byte c, float conf, float ras1, float ras2, float aa, float ab, float bb,
			float no) {
		probeSetName = psName;
		call = c;
		confidence = conf;
		RAS1 = ras1;
		RAS2 = ras2;
		aaCall = aa;
		abCall = ab;
		bbCall = bb;
		noCall = no;
	}

	public CHPGenotypeEntry(CHPGenotypeEntry entry) {
		probeSetName = entry.getProbeSetName();
		call = entry.getCall();
		confidence = entry.getConfidence();
		RAS1 = entry.getRAS1();
		RAS2 = entry.getRAS2();
		aaCall = entry.getAACall();
		abCall = entry.getABCall();
		bbCall = entry.getBBCall();
		noCall = entry.getNoCall();
	}

	public void clear() {
		probeSetName = null;
		call = 0;
		confidence = 0.0f;
		RAS1 = 0.0f;
		RAS2 = 0.0f;
		aaCall = 0.0f;
		abCall = 0.0f;
		bbCall = 0.0f;
		noCall = 0.0f;
	}

	public String getProbeSetName() {
		return probeSetName;
	}

	public byte getCall() {
		return call;
	}

	public float getConfidence() {
		return confidence;
	}

	public float getRAS1() {
		return RAS1;
	}

	public float getRAS2() {
		return RAS2;
	}

	public float getAACall() {
		return aaCall;
	}

	public float getABCall() {
		return abCall;
	}

	public float getBBCall() {
		return bbCall;
	}

	public float getNoCall() {
		return noCall;
	}

	public void setProbeSetName(String p) {
		probeSetName = p;
	}

	public void setCall(byte p) {
		call = p;
	}

	public void setConfidence(float p) {
		confidence = p;
	}

	public void setRAS1(float p) {
		RAS1 = p;
	}

	public void setRAS2(float p) {
		RAS2 = p;
	}

	public void setAACall(float p) {
		aaCall = p;
	}

	public void setABCall(float p) {
		abCall = p;
	}

	public void setBBCall(float p) {
		bbCall = p;
	}

	public void setNoCall(float p) {
		noCall = p;
	}
}
