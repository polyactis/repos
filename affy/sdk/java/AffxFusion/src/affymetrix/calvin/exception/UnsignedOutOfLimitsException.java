package affymetrix.calvin.exception;

public class UnsignedOutOfLimitsException extends CalvinException {
	public static final long serialVersionUID = 0;

	public UnsignedOutOfLimitsException(String description, int err) {
		super(description, err);
	}
}
