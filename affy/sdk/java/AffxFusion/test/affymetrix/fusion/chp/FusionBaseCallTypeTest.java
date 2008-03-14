/*
 * FusionBaseCallTypeTest.java
 * JUnit based test
 *
 * Created on December 2, 2005, 3:46 PM
 */

package affymetrix.fusion.chp;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class FusionBaseCallTypeTest extends TestCase {

	public FusionBaseCallTypeTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(FusionBaseCallTypeTest.class);

		return suite;
	}

	/**
	 * Test of setPosition method, of class affymetrix.fusion.chp.FusionBaseCallType.
	 */
	public void testProperties() {
		FusionBaseCallType b = new FusionBaseCallType();
		b.setCall((byte)'a');
		assertEquals(b.getCall(), 'a');
		b.setPosition(10);
		assertEquals(b.getPosition(), 10);

	}

}
