package affymetrix.calvin.data;

import affymetrix.portability.UByte;

public class GenotypeBufferEntry {

	private UByte call = UByte.ZERO;

	private float confidence = 0f;

	private float RAS1 = 0f;

	private float RAS2 = 0f;

	private float aaCall = 0f;

	private float abCall = 0f;

	private float bbCall = 0f;

	private float noCall = 0f;

	public GenotypeBufferEntry() {
	}

	public float getAaCall() {
		return aaCall;
	}

	public void setAaCall(float aaCall) {
		this.aaCall = aaCall;
	}

	public float getAbCall() {
		return abCall;
	}

	public void setAbCall(float abCall) {
		this.abCall = abCall;
	}

	public float getBbCall() {
		return bbCall;
	}

	public void setBbCall(float bbCall) {
		this.bbCall = bbCall;
	}

	public UByte getCall() {
		return call;
	}

	public void setCall(UByte call) {
		this.call = call;
	}

	public float getConfidence() {
		return confidence;
	}

	public void setConfidence(float confidence) {
		this.confidence = confidence;
	}

	public float getNoCall() {
		return noCall;
	}

	public void setNoCall(float noCall) {
		this.noCall = noCall;
	}

	public float getRAS1() {
		return RAS1;
	}

	public void setRAS1(float ras1) {
		RAS1 = ras1;
	}

	public float getRAS2() {
		return RAS2;
	}

	public void setRAS2(float ras2) {
		RAS2 = ras2;
	}
}
