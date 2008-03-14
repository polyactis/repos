package affymetrix.portability;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import junit.framework.TestCase;

public class UIntTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testZero() throws Exception {
		UInt uint = new UInt(23);
		assertTrue(uint.toLong() == 23);
		uint.zero();
		assertTrue(uint.toLong() == 0);
	}

	public void testSetLong() throws Exception {
		UInt uInt = new UInt();
		uInt.set(UInt.UINT_MAX - 1);
		assertTrue(uInt.toLong() == UInt.UINT_MAX - 1);
		ByteBuffer buf = ByteBuffer.allocate(8);
		buf.putLong(UInt.UINT_MAX - 1);
		byte[] b = buf.order(ByteOrder.BIG_ENDIAN).array();
		byte[] b2 = uInt.getBytes();
		for (int n = 4; n < 8; n++) {
			assertTrue(b[n] == b2[n - 4]);
		}
		for (long i = 0; i < UInt.UINT_MAX; i+=100) {
			uInt.zero();
			uInt.set(i);
			assertTrue(uInt.toLong() == i);
			buf = ByteBuffer.allocate(8);
			buf.putLong(i);
			b = buf.order(ByteOrder.BIG_ENDIAN).array();
			b2 = uInt.getBytes();
			for (int n = 4; n < 8; n++) {
				assertTrue(b[n] == b2[n - 4]);
			}
		}
	}

	public void testSetByteArray() throws Exception {
		ByteBuffer buf = ByteBuffer.allocate(4);
		buf.putInt(7000);
		byte[] b = buf.order(ByteOrder.BIG_ENDIAN).array();
		UInt uInt = new UInt();
		uInt.set(b);
		assertTrue(uInt.toLong() == 7000);
		byte[] b2 = uInt.getBytes();
		for (int i = 0; i < b.length; i++) {
			assertTrue(b[i] == b2[i]);
		}
	}

	public void testAddLong() throws Exception {
		UInt uInt = new UInt(39);
		uInt.add(3);
		assertTrue(uInt.toLong() == 42);
	}

	public void testSubtractLong() throws Exception {
		UInt uInt = new UInt(39);
		uInt.subtract(3);
		assertTrue(uInt.toLong() == 36);

		UInt uInt2 = new UInt(10);
		uInt.subtract(uInt2);
		assertTrue(uInt.toLong() == 26);
		assertTrue(uInt2.toLong() == 10);
	}

	public void testMultiplyLong() throws Exception {
		UInt uInt = new UInt(39);
		uInt.multiply(2);
		assertTrue(uInt.toLong() == 78);
	}

	public void testDivideLong() throws Exception {
		UInt uInt = new UInt(39);
		uInt.divide(2);
		assertTrue(uInt.toLong() == 19);
	}

	public void testGreaterThanLong() throws Exception {
		UInt uInt = new UInt(39);
		assertTrue(uInt.greaterThan(37));
		assertFalse(uInt.greaterThan(40));
	}

	public void testLessThanLong() throws Exception {
		UInt uInt = new UInt(39);
		assertTrue(uInt.lessThan(40));
		assertFalse(uInt.lessThan(37));
	}

	public void testGreaterThanOrEqualToLong() throws Exception {
		UInt uInt = new UInt(39);
		assertTrue(uInt.greaterThanOrEqualTo(37));
		assertTrue(uInt.greaterThanOrEqualTo(39));
		assertFalse(uInt.greaterThanOrEqualTo(40));
	}

	public void testLessThanOrEqualToLong() throws Exception {
		UInt uInt = new UInt(39);
		assertTrue(uInt.lessThanOrEqualTo(40));
		assertTrue(uInt.lessThanOrEqualTo(39));
		assertFalse(uInt.lessThanOrEqualTo(37));
	}

	public void testEqualsInt() throws Exception {
		UInt uInt = new UInt(39);
		assertTrue(uInt.equals(39));
		assertFalse(uInt.equals(37));
	}
}
