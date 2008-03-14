/*
 * ProbeSetInfoTest.java
 * JUnit based test
 *
 * Created on October 12, 2005, 2:22 PM
 */

package affymetrix.gcos.psi;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class ProbeSetInfoTest extends TestCase {

	public ProbeSetInfoTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(ProbeSetInfoTest.class);

		return suite;
	}

	/**
	 * Test of setProbeSetName method, of class affymetrix.gcos.psi.ProbeSetInfo.
	 */
	public void testProbeSetName() {
		ProbeSetInfo ps = new ProbeSetInfo();
		ps.setProbeSetName("ps");
		assertEquals(ps.getProbeSetName(), "ps");
	}

	/**
	 * Test of setNumberPairs method, of class affymetrix.gcos.psi.ProbeSetInfo.
	 */
	public void testNumberPairs() {
		ProbeSetInfo ps = new ProbeSetInfo();
		ps.setNumberPairs(10);
		assertEquals(ps.getNumberPairs(), 10);
	}

}
