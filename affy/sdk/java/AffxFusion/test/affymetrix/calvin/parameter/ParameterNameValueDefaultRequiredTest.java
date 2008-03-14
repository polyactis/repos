/*
 * ParameterNameValueDefaultRequiredTest.java
 * JUnit based test
 *
 * Created on October 23, 2005, 7:05 PM
 */

package affymetrix.calvin.parameter;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class ParameterNameValueDefaultRequiredTest extends TestCase {

	public ParameterNameValueDefaultRequiredTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(ParameterNameValueDefaultRequiredTest.class);

		return suite;
	}

	/**
	 * Test of ParameterValueTypeToString method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testParameterValueTypeToString() {
	}

	/**
	 * Test of ParameterValueTypeFromString method, of class
	 * affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testParameterValueTypeFromString() {
	}

	/**
	 * Test of getControlMultiValues method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testGetControlMultiValues() {
	}

	/**
	 * Test of getDefaultValueType method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testgetDefaultValueType() {
	}

	/**
	 * Test of ControlledVocabulary method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testControlledVocabulary() {
	}

	/**
	 * Test of getDefaultMIMEValue method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testGetDefaultMIMEValue() {
	}

	/**
	 * Test of getRequiredFlag method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testRequiredFlag() {
		ParameterNameValueDefaultRequired p = new ParameterNameValueDefaultRequired();
		assertEquals(p.getRequiredFlag(), false);
		p.setRequiredFlag(true);
		assertEquals(p.getRequiredFlag(), true);
	}

	/**
	 * Test of getHasDefault method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testHasDefault() {
		ParameterNameValueDefaultRequired p = new ParameterNameValueDefaultRequired();
		p.setName("p");
		p.setValueInt32(123);
		assertEquals(p.getHasDefault(), false);
		p.setDefaultValueInt32(321);
		assertEquals(p.getHasDefault(), true);
	}

	/**
	 * Test of defaultToString method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testDefaultToString() {
		// ParameterNameValueDefaultRequired p = new ParameterNameValueDefaultRequired();
	}

	/**
	 * Test of getDefaultValueInt8 method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testValueInt8() {
		ParameterNameValueDefaultRequired p = new ParameterNameValueDefaultRequired();
		p.setName("p");
		p.setDefaultValueInt8((byte)12);
		assertEquals(p.getDefaultValueInt8(), (byte)12);
	}

	/**
	 * Test of getDefaultValueUInt8 method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testValueUInt8() {
		ParameterNameValueDefaultRequired p = new ParameterNameValueDefaultRequired();
		p.setName("p");
		p.setDefaultValueUInt8((byte)12);
		assertEquals(p.getDefaultValueUInt8(), (byte)12);
	}

	/**
	 * Test of getDefaultValueInt16 method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testValueInt16() {
		ParameterNameValueDefaultRequired p = new ParameterNameValueDefaultRequired();
		p.setName("p");
		p.setDefaultValueInt16((short)12);
		assertEquals(p.getDefaultValueInt16(), (short)12);
	}

	/**
	 * Test of getDefaultValueUInt16 method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testValueUInt16() {
		ParameterNameValueDefaultRequired p = new ParameterNameValueDefaultRequired();
		p.setName("p");
		p.setDefaultValueUInt16((short)12);
		assertEquals(p.getDefaultValueUInt16(), (short)12);
	}

	/**
	 * Test of getDefaultValueInt32 method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testValueInt32() {
		ParameterNameValueDefaultRequired p = new ParameterNameValueDefaultRequired();
		p.setName("p");
		p.setDefaultValueInt32(1);
		assertEquals(1, p.getDefaultValueInt32());
		p.setDefaultValueInt32(2);
		assertEquals(2, p.getDefaultValueInt32());
		p.setDefaultValueInt32(127);
		assertEquals(127, p.getDefaultValueInt32());
		p.setDefaultValueInt32(128);
		assertEquals(128, p.getDefaultValueInt32());
		p.setDefaultValueInt32(255);
		assertEquals(255, p.getDefaultValueInt32());
		p.setDefaultValueInt32(256);
		assertEquals(256, p.getDefaultValueInt32());
		p.setDefaultValueInt32(1024);
		assertEquals(1024, p.getDefaultValueInt32());
		p.setDefaultValueInt32(10240);
		assertEquals(10240, p.getDefaultValueInt32());
		p.setDefaultValueInt32(102450);
		assertEquals(102450, p.getDefaultValueInt32());
		p.setDefaultValueInt32(Integer.MAX_VALUE - 1);
		assertEquals(p.getDefaultValueInt32(), Integer.MAX_VALUE - 1);
	}

	/**
	 * Test of getDefaultValueUInt32 method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testValueUInt32() {
		ParameterNameValueDefaultRequired p = new ParameterNameValueDefaultRequired();
		p.setName("p");
		p.setDefaultValueUInt32(12);
		assertEquals(p.getDefaultValueUInt32(), 12);
	}

	/**
	 * Test of getDefaultValueFloat method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testValueFloat() {
		ParameterNameValueDefaultRequired p = new ParameterNameValueDefaultRequired();
		p.setName("p");
		p.setDefaultValueFloat(12121212.123f);
		assertEquals(p.getDefaultValueFloat(), 12121212.123f, 0.001f);
		p.setDefaultValueFloat(121.123f);
		assertEquals(p.getDefaultValueFloat(), 121.123f, 0.001f);
		p.setDefaultValueFloat(122.123f);
		assertEquals(p.getDefaultValueFloat(), 122.123f, 0.001f);
		p.setDefaultValueFloat(1233.123f);
		assertEquals(p.getDefaultValueFloat(), 1233.123f, 0.001f);
		p.setDefaultValueFloat(12334.123f);
		assertEquals(p.getDefaultValueFloat(), 12334.123f, 0.001f);
	}

	/**
	 * Test of getDefaultValueText method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testValueText() {
		ParameterNameValueDefaultRequired p = new ParameterNameValueDefaultRequired();
		p.setName("p");
		p.setDefaultValueText("text");
		assertEquals(p.getDefaultValueText(), "text");
	}

	/**
	 * Test of getDefaultValueAscii method, of class affymetrix.calvin.parameter.ParameterNameValueDefaultRequired.
	 */
	public void testValueAscii() {
		ParameterNameValueDefaultRequired p = new ParameterNameValueDefaultRequired();
		p.setName("p");
		p.setDefaultValueAscii("text");
		assertEquals(p.getDefaultValueAscii(), "text");
	}

}
