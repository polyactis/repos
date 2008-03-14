/*
 * BaseCallTypeTest.java
 * JUnit based test
 *
 * Created on October 13, 2005, 1:53 PM
 */

package affymetrix.gcos.chp;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class BaseCallTypeTest extends TestCase {

	public BaseCallTypeTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(BaseCallTypeTest.class);

		return suite;
	}

	/**
	 * Test of properties.
	 */
	public void testProperties() {
		BaseCallType t = new BaseCallType();
		assertEquals(t.getPosition(), 0);
		assertEquals(t.getCall(), ' ');

		t.setPosition(10);
		t.setCall('a');

		assertEquals(t.getPosition(), 10);
		assertEquals(t.getCall(), 'a');

	}

	/**
	 * Test of the constructors.
	 */
	public void testConstructor() {
		BaseCallType t1 = new BaseCallType();

		assertEquals(t1.getPosition(), 0);
		assertEquals(t1.getCall(), ' ');

		t1.setPosition(10);
		t1.setCall('a');

		BaseCallType t2 = new BaseCallType(t1);
		t1.setPosition(0);
		t1.setCall(' ');

		assertEquals(t2.getPosition(), 10);
		assertEquals(t2.getCall(), 'a');
	}
}
