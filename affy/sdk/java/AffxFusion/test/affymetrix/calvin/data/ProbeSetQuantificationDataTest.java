/*
 * ProbeSetQuantificationDataTest.java
 * JUnit based test
 *
 * Created on December 23, 2005, 6:36 PM
 */

package affymetrix.calvin.data;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class ProbeSetQuantificationDataTest extends TestCase {

	public ProbeSetQuantificationDataTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(ProbeSetQuantificationDataTest.class);

		return suite;
	}

	/**
	 * Test of getName method, of class affymetrix.calvin.data.ProbeSetQuantificationData.
	 */
	public void testProperties() {
		ProbeSetQuantificationData d = new ProbeSetQuantificationData();
		d.setQuantification(10.0f);
		assertEquals(d.getQuantification(), 10.0f, 0.0001f);
		d.setName("abc");
		assertEquals(d.getName(), "abc");
		d.setId(10);
		assertEquals(d.getId(), 10);
	}

}
