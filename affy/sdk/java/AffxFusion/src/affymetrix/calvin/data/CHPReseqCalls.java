package affymetrix.calvin.data;

public class CHPReseqCalls {

	/** The force call was made due to no signal threshold. */
	public static final byte CC_NO_SIGNAL_THR_FORCE_CALL = 'N';

	/** The force call was made due to weak signal threshold. */
	public static final byte CC_WEAK_SIGNAL_THR_FORCE_CALL = 'W';

	/** The force call was made due to saturation level. */
	public static final byte CC_SATURATION_LEVEL_FORCE_CALL = 'S';

	/** The force call was made due to quality score threshold. */
	public static final byte CC_QUALITY_SCORE_THR_FORCE_CALL = 'Q';

	/** The force call was made due to failed both trace and sequence profiles. */
	public static final byte CC_TRACE_AND_SEQUENCE_PROFILES_FORCE_CALL = 'F';

	/** The force call was made due to base reliability threshold. */
	public static final byte CC_RELIABILITY_THR_FORCE_CALL = 'B';

	private CHPReseqCalls() {
	}

}
