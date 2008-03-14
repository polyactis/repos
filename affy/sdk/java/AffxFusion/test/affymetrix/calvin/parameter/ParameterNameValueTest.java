/*
 * ParameterNameValueTest.java
 * JUnit based test
 *
 * Created on October 23, 2005, 7:05 PM
 */

package affymetrix.calvin.parameter;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.portability.UByte;
import affymetrix.portability.UInt;
import affymetrix.portability.UShort;

/**
 * 
 * @author ljevon
 */
public class ParameterNameValueTest extends TestCase {

	public ParameterNameValueTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(ParameterNameValueTest.class);

		return suite;
	}

	/**
	 * Test of equals method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testEquals() throws Exception {
		ParameterNameValue p = new ParameterNameValue();
		p.setName("p");
		p.setValueAscii("val");
		assertEquals(p.equals("p"), true);
		assertEquals(p.equals("k"), false);
		ParameterNameValue p2 = new ParameterNameValue();
		p2.setName("p");
		assertEquals(p.equals(p2), true);
		p2.setName("l");
		assertEquals(p.equals(p2), false);
	}

	/**
	 * Test of getName method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testName() {

	}

	/**
	 * Test of getParameterType method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testGetParameterType() {

	}

	/**
	 * Test of getValueInt8 method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testValueInt8() {
		ParameterNameValue p = new ParameterNameValue();
		p.setName("p");
		p.setValueInt8((byte)12);
		assertEquals(p.getValueInt8(), (byte)12);
	}

	/**
	 * Test of getValueUInt8 method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testValueUInt8() throws Exception {
		ParameterNameValue p = new ParameterNameValue();
		p.setName("p");
		p.setValueUInt8(new UByte((short)12));
		assertEquals(p.getValueUInt8().toShort(), (short)12);
	}

	/**
	 * Test of getValueInt16 method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testValueInt16() {
		ParameterNameValue p = new ParameterNameValue();
		p.setName("p");
		p.setValueInt16((short)12);
		assertEquals(p.getValueInt16(), (short)12);
	}

	/**
	 * Test of getValueUInt16 method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testValueUInt16() throws Exception {
		ParameterNameValue p = new ParameterNameValue();
		p.setName("p");
		p.setValueUInt16(new UShort(12));
		assertEquals(p.getValueUInt16().toInt(), 12);
	}

	/**
	 * Test of getValueInt32 method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testValueInt32() {
		ParameterNameValue p = new ParameterNameValue();
		p.setName("p");
		p.setValueInt32(1);
		assertEquals(1, p.getValueInt32());
		p.setValueInt32(2);
		assertEquals(2, p.getValueInt32());
		p.setValueInt32(127);
		assertEquals(127, p.getValueInt32());
		p.setValueInt32(128);
		assertEquals(128, p.getValueInt32());
		p.setValueInt32(255);
		assertEquals(255, p.getValueInt32());
		p.setValueInt32(256);
		assertEquals(256, p.getValueInt32());
		p.setValueInt32(1024);
		assertEquals(1024, p.getValueInt32());
		p.setValueInt32(10240);
		assertEquals(10240, p.getValueInt32());
		p.setValueInt32(102450);
		assertEquals(102450, p.getValueInt32());
		p.setValueInt32(Integer.MAX_VALUE - 1);
		assertEquals(p.getValueInt32(), Integer.MAX_VALUE - 1);
	}

	/**
	 * Test of getValueUInt32 method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testValueUInt32() throws Exception {
		ParameterNameValue p = new ParameterNameValue();
		p.setName("p");
		p.setValueUInt32(new UInt(12));
		assertEquals(p.getValueUInt32().toLong(), 12);
	}

	/**
	 * Test of getValueFloat method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testValueFloat() {
		ParameterNameValue p = new ParameterNameValue();
		p.setName("p");
		p.setValueFloat(12121212.123f);
		assertEquals(p.getValueFloat(), 12121212.123f, 0.001f);
		p.setValueFloat(121.123f);
		assertEquals(p.getValueFloat(), 121.123f, 0.001f);
		p.setValueFloat(122.123f);
		assertEquals(p.getValueFloat(), 122.123f, 0.001f);
		p.setValueFloat(1233.123f);
		assertEquals(p.getValueFloat(), 1233.123f, 0.001f);
		p.setValueFloat(12334.123f);
		assertEquals(p.getValueFloat(), 12334.123f, 0.001f);
	}

	/**
	 * Test of getValueText method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testValueText() {
		ParameterNameValue p = new ParameterNameValue();
		p.setName("p");
		p.setValueText("text");
		assertEquals(p.getValueText(), "text");
	}

	/**
	 * Test of getValueAscii method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testValueAscii() throws Exception {
		ParameterNameValue p = new ParameterNameValue();
		p.setName("p");
		p.setValueAscii("text");
		assertEquals(p.getValueAscii(), "text");
	}

	/**
	 * Test of getMIMEType method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testMIMEType() {
	}

	/**
	 * Test of getMIMEValue method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testMIMEValue() {

	}

	/**
	 * Test of toString method, of class affymetrix.calvin.parameter.ParameterNameValue.
	 */
	public void testToString() {

	}

}
