package affymetrix.portability;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import junit.framework.TestCase;

public class UByteTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testZero() throws Exception {
		UByte uByte = new UByte(((short)23));
		assertTrue(uByte.toShort() == 23);
		uByte.zero();
		assertTrue(uByte.toShort() == 0);
	}

	public void testSetShort() throws Exception {
		UByte uByte = new UByte();
		for (int i = 0; i < UByte.UBYTE_MAX; i++) {
			uByte.zero();
			uByte.set(((short)i));
			assertTrue(uByte.toShort() == i);
			ByteBuffer buf = ByteBuffer.allocate(1);
			buf.put(((byte)i));
			byte[] b = buf.order(ByteOrder.BIG_ENDIAN).array();
			byte[] b2 = uByte.getBytes();
			for (int n = 0; n < b.length; n++) {
				assertTrue(b[n] == b2[n]);
			}
		}
	}

	public void testSetByteArray() throws Exception {
		ByteBuffer buf = ByteBuffer.allocate(1);
		buf.put(((byte)44));
		byte[] b = buf.order(ByteOrder.BIG_ENDIAN).array();
		UByte uByte = new UByte();
		uByte.set(b);
		assertTrue(uByte.toShort() == 44);
		byte[] b2 = uByte.getBytes();
		for (int i = 0; i < b.length; i++) {
			assertTrue(b[i] == b2[i]);
		}
	}

	public void testAddLong() throws Exception {
		UByte uByte = new UByte(((short)39));
		uByte.add(((short)3));
		assertTrue(uByte.toShort() == 42);
	}

	public void testSubtractLong() throws Exception {
		UByte uByte = new UByte(((short)39));
		uByte.subtract(((short)3));
		assertTrue(uByte.toShort() == 36);
	}

	public void testMultiplyLong() throws Exception {
		UByte uByte = new UByte(((short)39));
		uByte.multiply(((short)2));
		assertTrue(uByte.toShort() == 78);
	}

	public void testDivideLong() throws Exception {
		UByte uByte = new UByte(((short)39));
		uByte.divide(((short)2));
		assertTrue(uByte.toShort() == 19);
	}

	public void testGreaterThanLong() throws Exception {
		UByte uByte = new UByte(((short)39));
		assertTrue(uByte.greaterThan(((short)37)));
		assertFalse(uByte.greaterThan(((short)40)));
	}

	public void testLessThanLong() throws Exception {
		UByte uByte = new UByte(((short)39));
		assertTrue(uByte.lessThan(((short)40)));
		assertFalse(uByte.lessThan(((short)37)));
	}

	public void testGreaterThanOrEqualToLong() throws Exception {
		UByte uByte = new UByte(((short)39));
		assertTrue(uByte.greaterThanOrEqualTo(((short)37)));
		assertTrue(uByte.greaterThanOrEqualTo(((short)39)));
		assertFalse(uByte.greaterThanOrEqualTo(((short)40)));
	}

	public void testLessThanOrEqualToLong() throws Exception {
		UByte uByte = new UByte(((short)39));
		assertTrue(uByte.lessThanOrEqualTo(((short)40)));
		assertTrue(uByte.lessThanOrEqualTo(((short)39)));
		assertFalse(uByte.lessThanOrEqualTo(((short)37)));
	}

	public void testEqualsInt() throws Exception {
		UByte uByte = new UByte(((short)39));
		assertTrue(uByte.equals(((short)39)));
		assertFalse(uByte.equals(((short)37)));
	}

}
