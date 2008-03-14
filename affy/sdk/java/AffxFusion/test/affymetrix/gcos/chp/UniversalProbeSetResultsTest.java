/*
 * UniversalProbeSetResultsTest.java
 * JUnit based test
 *
 * Created on October 13, 2005, 9:19 AM
 */

package affymetrix.gcos.chp;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class UniversalProbeSetResultsTest extends TestCase {

	public UniversalProbeSetResultsTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(UniversalProbeSetResultsTest.class);

		return suite;
	}

	/**
	 * Test of getBackground method, of class affymetrix.gcos.chp.UniversalProbeSetResults.
	 */
	public void testBackground() {
		UniversalProbeSetResults u = new UniversalProbeSetResults();
		assertEquals(u.getBackground(), 0.0f, 0.0000001);
		u.setBackground(1.0f);
		assertEquals(u.getBackground(), 1.0f, 0.0000001);

	}

}
