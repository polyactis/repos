/*
 * ForceCallTypeTest.java
 * JUnit based test
 *
 * Created on October 13, 2005, 1:45 PM
 */

package affymetrix.gcos.chp;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class ForceCallTypeTest extends TestCase {

	public ForceCallTypeTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(ForceCallTypeTest.class);

		return suite;
	}

	/**
	 * Test of getPosition method, of class affymetrix.gcos.chp.ForceCallType.
	 */
	public void testProperties() throws Exception {
		ForceCallType t = new ForceCallType();
		assertEquals(t.getPosition(), 0);
		assertEquals(t.getCall(), (byte)' ');
		assertEquals(t.getReason(), (byte)' ');
		t.setPosition(10);
		t.setCall((byte)'a');
		t.setReason(ForceCallType.NO_SIGNAL_THR_FORCE_CALL);
		assertEquals(t.getPosition(), 10);
		assertEquals(t.getCall(), (byte)'a');
		assertEquals(t.getReason(), ForceCallType.NO_SIGNAL_THR_FORCE_CALL);

	}

	/**
	 * Test of the constructors.
	 */
	public void testConstructor() {
		ForceCallType t1 = new ForceCallType();
		assertEquals(t1.getPosition(), 0);
		assertEquals(t1.getCall(), (byte)' ');
		assertEquals(t1.getReason(), (byte)' ');
		t1.setPosition(10);
		t1.setCall((byte)'a');
		t1.setReason(ForceCallType.NO_SIGNAL_THR_FORCE_CALL);
		ForceCallType t2 = new ForceCallType(t1);
		t1.setPosition(0);
		t1.setReason(Byte.MIN_VALUE);
		t1.setCall((byte)' ');
		assertEquals(t2.getPosition(), 10);
		assertEquals(t2.getCall(), (byte)'a');
		assertEquals(t2.getReason(), ForceCallType.NO_SIGNAL_THR_FORCE_CALL);
	}
}
