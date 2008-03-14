package affymetrix.portability;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import junit.framework.TestCase;

public class UShortTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testZero() throws Exception {
		UShort uint = new UShort(23);
		assertTrue(uint.toInt() == 23);
		uint.zero();
		assertTrue(uint.toInt() == 0);
	}

	public void testSetInt() throws Exception {
		UShort uShort = new UShort();
		for (int i = 0; i < UShort.USHORT_MAX; i++) {
			uShort.zero();
			uShort.set(i);
			assertTrue(uShort.toInt() == i);
			ByteBuffer buf = ByteBuffer.allocate(2);
			buf.putShort(((short)i));
			byte[] b = buf.order(ByteOrder.BIG_ENDIAN).array();
			byte[] b2 = uShort.getBytes();
			for (int n = 0; n < b.length; n++) {
				assertTrue(b[n] == b2[n]);
			}
		}
	}

	public void testSetByteArray() throws Exception {
		ByteBuffer buf = ByteBuffer.allocate(2);
		buf.putShort(((short)44));
		byte[] b = buf.order(ByteOrder.BIG_ENDIAN).array();
		UShort uShort = new UShort();
		uShort.set(b);
		assertTrue(uShort.toInt() == 44);
		byte[] b2 = uShort.getBytes();
		for (int i = 0; i < b.length; i++) {
			assertTrue(b[i] == b2[i]);
		}
	}

	public void testAddLong() throws Exception {
		UShort uShort = new UShort(39);
		uShort.add(3);
		assertTrue(uShort.toInt() == 42);
	}

	public void testSubtractLong() throws Exception {
		UShort uShort = new UShort(39);
		uShort.subtract(3);
		assertTrue(uShort.toInt() == 36);
	}

	public void testMultiplyLong() throws Exception {
		UShort uShort = new UShort(39);
		uShort.multiply(2);
		assertTrue(uShort.toInt() == 78);
	}

	public void testDivideLong() throws Exception {
		UShort uShort = new UShort(39);
		uShort.divide(2);
		assertTrue(uShort.toInt() == 19);
	}

	public void testGreaterThanLong() throws Exception {
		UShort uShort = new UShort(39);
		assertTrue(uShort.greaterThan(37));
		assertFalse(uShort.greaterThan(40));
	}

	public void testLessThanLong() throws Exception {
		UShort uShort = new UShort(39);
		assertTrue(uShort.lessThan(40));
		assertFalse(uShort.lessThan(37));
	}

	public void testGreaterThanOrEqualToLong() throws Exception {
		UShort uShort = new UShort(39);
		assertTrue(uShort.greaterThanOrEqualTo(37));
		assertTrue(uShort.greaterThanOrEqualTo(39));
		assertFalse(uShort.greaterThanOrEqualTo(40));
	}

	public void testLessThanOrEqualToLong() throws Exception {
		UShort uShort = new UShort(39);
		assertTrue(uShort.lessThanOrEqualTo(40));
		assertTrue(uShort.lessThanOrEqualTo(39));
		assertFalse(uShort.lessThanOrEqualTo(37));
	}

	public void testEqualsInt() throws Exception {
		UShort uShort = new UShort(39);
		assertTrue(uShort.equals(39));
		assertFalse(uShort.equals(37));
	}
}
