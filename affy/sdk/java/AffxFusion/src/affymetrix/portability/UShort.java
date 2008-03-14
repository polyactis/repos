package affymetrix.portability;

import java.util.List;

import affymetrix.calvin.exception.UnsignedOutOfLimitsException;

public class UShort extends UnsignedBase {

	public static final int USHORT_MAX = 65536;

	public static final UShort ZERO = new UShort();

	private int uShort = 0;

	public UShort() {
	}

	public UShort(int l) throws UnsignedOutOfLimitsException {
		this();
		set(l);
	}

	public UShort(byte[] b) throws UnsignedOutOfLimitsException {
		this();
		set(b);
	}

	@Override
	public void zero() {
		uShort = 0;
	}

	public void set(int l) throws UnsignedOutOfLimitsException {
		if (inLimits(l, USHORT_MAX)) {
			uShort = l;
		}
	}

	public void set(byte[] b) throws UnsignedOutOfLimitsException {
		if (inLimits(b, DataSizes.SHORT_SIZE)) {
			uShort = (int)unsignedBytesToLong(b);
		}
	}

	public int add(UShort rhs) throws UnsignedOutOfLimitsException {
		return add(rhs.toInt());
	}

	public int add(int rhs) throws UnsignedOutOfLimitsException {
		int i = uShort + rhs;
		if (inLimits(i, USHORT_MAX)) {
			uShort = i;
		}
		return uShort;
	}

	public int subtract(UShort rhs) throws UnsignedOutOfLimitsException {
		return subtract(rhs.toInt());
	}

	public int subtract(int rhs) throws UnsignedOutOfLimitsException {
		int i = uShort - rhs;
		if (inLimits(i, USHORT_MAX)) {
			uShort = i;
		}
		return uShort;
	}

	public int multiply(UShort rhs) throws UnsignedOutOfLimitsException {
		return multiply(rhs.toInt());
	}

	public int multiply(int rhs) throws UnsignedOutOfLimitsException {
		int i = uShort * rhs;
		if (inLimits(i, USHORT_MAX)) {
			uShort = i;
		}
		return uShort;
	}

	public int divide(UShort rhs) throws UnsignedOutOfLimitsException {
		return divide(rhs.toInt());
	}

	public int divide(int rhs) throws UnsignedOutOfLimitsException {
		int i = uShort / rhs;
		if (inLimits(i, USHORT_MAX)) {
			uShort = i;
		}
		return uShort;
	}

	public boolean greaterThan(UShort rhs) {
		return greaterThan(rhs.toInt());
	}

	public boolean greaterThan(int rhs) {
		return uShort > rhs;
	}

	public boolean lessThan(UShort rhs) {
		return lessThan(rhs.toInt());
	}

	public boolean lessThan(int rhs) {
		return uShort < rhs;
	}

	public boolean greaterThanOrEqualTo(UShort rhs) {
		return greaterThanOrEqualTo(rhs.toInt());
	}

	public boolean greaterThanOrEqualTo(int rhs) {
		return uShort >= rhs;
	}

	public boolean lessThanOrEqualTo(UShort rhs) {
		return lessThanOrEqualTo(rhs.toInt());
	}

	public boolean lessThanOrEqualTo(int rhs) {
		return uShort <= rhs;
	}

	public boolean equals(UShort rhs) {
		return equals(rhs.toInt());
	}

	public boolean equals(int rhs) {
		return uShort == rhs;
	}

	public int toInt() {
		return uShort;
	}

	@Override
	public byte[] getBytes() {
		return toUnsignedByteArray(uShort);
	}

	public static byte[] toUnsignedByteArray(int value) {
		return toUnsignedByteArray(value, DataSizes.SHORT_SIZE);
	}

	public static void addToUnsignedByteList(int value, List<Byte> buf) {
		addBytes(toUnsignedByteArray(value), buf);
	}
}
