package affymetrix.calvin.utils;

public class IOUtils {

	public static final String EMPTY = "";

	public static final String UNICODE_CHARSET = "UTF-16BE";

	public static final String ASCII_CHARSET = "US-ASCII";

	private IOUtils() {
	}

	public static boolean isNullOrEmpty(String s) {
		if ((s == null) || (s.length() == 0)) {
			return true;
		}
		return false;
	}

	public static String getVal(String val) {
		return getVal(val, true);
	}

	public static String getVal(String val, boolean trim) {
		return getVal(val, EMPTY, trim);
	}

	public static String getVal(String val, String defaultVal, boolean trim) {
		if (isNullOrEmpty(val)) {
			return defaultVal;
		}
		if (trim) {
			return val.trim();
		}
		return val;
	}

}
