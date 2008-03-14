package affymetrix.calvin.data;

import affymetrix.portability.UInt;

public class ProbeSetMultiDataCytoRegionData extends ProbeSetMultiDataBase {

	/** A no call for the cyto result. */
	public static final byte CYTO_NO_CALL = 0;

	/** An absent call for the cyto result. */
	public static final byte CYTO_ABSENT_CALL = 1;

	/** A present call for the cyto result. */
	public static final byte CYTO_PRESENT_CALL = 2;

	/** The name of the region. */
	private String name = null;

	/** The chromosome value. */
	private byte chr = 0;

	/** The physical start position. */
	private UInt startPosition = null;

	/** The physical stop position. */
	private UInt stopPosition = null;

	/** The call for the region. */
	private byte call = 0;

	/** The confidence score. */
	private float confidence = 0f;

	public byte getCall() {
		return call;
	}

	public void setCall(byte call) {
		this.call = call;
	}

	public byte getChr() {
		return chr;
	}

	public void setChr(byte chr) {
		this.chr = chr;
	}

	public float getConfidence() {
		return confidence;
	}

	public void setConfidence(float confidence) {
		this.confidence = confidence;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public UInt getStartPosition() {
		return startPosition;
	}

	public void setStartPosition(UInt startPosition) {
		this.startPosition = startPosition;
	}

	public UInt getStopPosition() {
		return stopPosition;
	}

	public void setStopPosition(UInt stopPosition) {
		this.stopPosition = stopPosition;
	}
}
