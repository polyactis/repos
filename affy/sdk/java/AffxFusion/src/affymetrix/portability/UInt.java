package affymetrix.portability;

import java.util.List;

import affymetrix.calvin.exception.UnsignedOutOfLimitsException;

public class UInt extends UnsignedBase {

	public static final long UINT_MAX = 4294967296L;

	public static final UInt ZERO = new UInt();

	private long uInt = 0;

	public UInt() {
	}

	public UInt(long l) throws UnsignedOutOfLimitsException {
		this();
		set(l);
	}

	public UInt(byte[] b) throws UnsignedOutOfLimitsException {
		this();
		set(b);
	}

	public UInt(UInt b) throws UnsignedOutOfLimitsException {
		this();
		set(b.toLong());
	}

	@Override
	public void zero() {
		uInt = 0;
	}

	public void set(long l) throws UnsignedOutOfLimitsException {
		if (inLimits(l, UINT_MAX)) {
			uInt = l;
		}
	}

	public void set(byte[] b) throws UnsignedOutOfLimitsException {
		if (inLimits(b, DataSizes.INT_SIZE)) {
			uInt = unsignedBytesToLong(b);
		}
	}

	public long add(UInt rhs) throws UnsignedOutOfLimitsException {
		return add(rhs.toLong());
	}

	public long add(long rhs) throws UnsignedOutOfLimitsException {
		long l = uInt + rhs;
		if (inLimits(l, UINT_MAX)) {
			uInt = l;
		}
		return uInt;
	}

	public long subtract(UInt rhs) throws UnsignedOutOfLimitsException {
		return subtract(rhs.toLong());
	}

	public long subtract(long rhs) throws UnsignedOutOfLimitsException {
		long l = uInt - rhs;
		if (inLimits(l, UINT_MAX)) {
			uInt = l;
		}
		return uInt;
	}

	public long multiply(UInt rhs) throws UnsignedOutOfLimitsException {
		return multiply(rhs.toLong());
	}

	public long multiply(long rhs) throws UnsignedOutOfLimitsException {
		long l = uInt * rhs;
		if (inLimits(l, UINT_MAX)) {
			uInt = l;
		}
		return uInt;
	}

	public long divide(UInt rhs) throws UnsignedOutOfLimitsException {
		return divide(rhs.toLong());
	}

	public long divide(long rhs) throws UnsignedOutOfLimitsException {
		long l = uInt / rhs;
		if (inLimits(l, UINT_MAX)) {
			uInt = l;
		}
		return uInt;
	}

	public boolean greaterThan(UInt rhs) {
		return greaterThan(rhs.toLong());
	}

	public boolean greaterThan(long rhs) {
		return uInt > rhs;
	}

	public boolean lessThan(UInt rhs) {
		return lessThan(rhs.toLong());
	}

	public boolean lessThan(long rhs) {
		return uInt < rhs;
	}

	public boolean greaterThanOrEqualTo(UInt rhs) {
		return greaterThanOrEqualTo(rhs.toLong());
	}

	public boolean greaterThanOrEqualTo(long rhs) {
		return uInt >= rhs;
	}

	public boolean lessThanOrEqualTo(UInt rhs) {
		return lessThanOrEqualTo(rhs.toLong());
	}

	public boolean lessThanOrEqualTo(long rhs) {
		return uInt <= rhs;
	}

	public boolean equals(UInt rhs) {
		return equals(rhs.toLong());
	}

	public boolean equals(long rhs) {
		return uInt == rhs;
	}

	public long toLong() {
		return uInt;
	}

	@Override
	public byte[] getBytes() {
		return toUnsignedByteArray(uInt);
	}

	public static byte[] toUnsignedByteArray(long value) {
		return toUnsignedByteArray(value, DataSizes.INT_SIZE);
	}

	public static void addToUnsignedByteList(long value, List<Byte> buf) {
		addBytes(toUnsignedByteArray(value), buf);
	}
}
