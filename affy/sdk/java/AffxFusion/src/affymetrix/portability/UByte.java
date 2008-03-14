package affymetrix.portability;

import java.util.List;

import affymetrix.calvin.exception.UnsignedOutOfLimitsException;

public class UByte extends UnsignedBase {

	public static final int UBYTE_MAX = 256;

	public static final UByte ZERO = new UByte();

	private short uByte = 0;

	public UByte() {
	}

	public UByte(short l) throws UnsignedOutOfLimitsException {
		this();
		set(l);
	}

	public UByte(byte[] b) throws UnsignedOutOfLimitsException {
		this();
		set(b);
	}

	public void set(short l) throws UnsignedOutOfLimitsException {
		if (inLimits(l, UBYTE_MAX)) {
			uByte = l;
		}
	}

	@Override
	public void zero() {
		uByte = 0;
	}

	public void set(byte[] b) throws UnsignedOutOfLimitsException {
		if (inLimits(b, DataSizes.CHAR_SIZE)) {
			uByte = (short)unsignedBytesToLong(b);
		}
	}

	public short add(UByte rhs) throws UnsignedOutOfLimitsException {
		return add(rhs.toShort());
	}

	public short add(short rhs) throws UnsignedOutOfLimitsException {
		int i = uByte + rhs;
		if (inLimits(i, UBYTE_MAX)) {
			uByte = (short)i;
		}
		return uByte;
	}

	public short subtract(UByte rhs) throws UnsignedOutOfLimitsException {
		return subtract(rhs.toShort());
	}

	public short subtract(short rhs) throws UnsignedOutOfLimitsException {
		int i = uByte - rhs;
		if (inLimits(i, UBYTE_MAX)) {
			uByte = (short)i;
		}
		return uByte;
	}

	public short multiply(UByte rhs) throws UnsignedOutOfLimitsException {
		return multiply(rhs.toShort());
	}

	public short multiply(short rhs) throws UnsignedOutOfLimitsException {
		int i = uByte * rhs;
		if (inLimits(i, UBYTE_MAX)) {
			uByte = (short)i;
		}
		return uByte;
	}

	public short divide(UByte rhs) throws UnsignedOutOfLimitsException {
		return divide(rhs.toShort());
	}

	public short divide(short rhs) throws UnsignedOutOfLimitsException {
		int i = uByte / rhs;
		if (inLimits(i, UBYTE_MAX)) {
			uByte = (short)i;
		}
		return uByte;
	}

	public boolean greaterThan(UByte rhs) {
		return greaterThan(rhs.toShort());
	}

	public boolean greaterThan(short rhs) {
		return uByte > rhs;
	}

	public boolean lessThan(UByte rhs) {
		return lessThan(rhs.toShort());
	}

	public boolean lessThan(short rhs) {
		return uByte < rhs;
	}

	public boolean greaterThanOrEqualTo(UByte rhs) {
		return greaterThanOrEqualTo(rhs.toShort());
	}

	public boolean greaterThanOrEqualTo(short rhs) {
		return uByte >= rhs;
	}

	public boolean lessThanOrEqualTo(UByte rhs) {
		return lessThanOrEqualTo(rhs.toShort());
	}

	public boolean lessThanOrEqualTo(short rhs) {
		return uByte <= rhs;
	}

	public boolean equals(UByte rhs) {
		return equals(rhs.toShort());
	}

	public boolean equals(short rhs) {
		return uByte == rhs;
	}

	public short toShort() {
		return uByte;
	}

	@Override
	public byte[] getBytes() {
		return toUnsignedByteArray(uByte);
	}

	public static byte[] toUnsignedByteArray(short value) {
		return toUnsignedByteArray(value, DataSizes.CHAR_SIZE);
	}

	public static void addToUnsignedByteList(short value, List<Byte> buf) {
		addBytes(toUnsignedByteArray(value), buf);
	}
}
