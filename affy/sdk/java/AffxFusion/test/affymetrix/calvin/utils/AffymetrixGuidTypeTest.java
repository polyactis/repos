/*
 * AffymetrixGuidTypeTest.java
 * JUnit based test
 *
 * Created on October 23, 2005, 7:05 PM
 */

package affymetrix.calvin.utils;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class AffymetrixGuidTypeTest extends TestCase {

	public AffymetrixGuidTypeTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(AffymetrixGuidTypeTest.class);

		return suite;
	}

	/**
	 * Test of getGuid method, of class affymetrix.calvin.utils.AffymetrixGuidType.
	 */
	public void testGuid() throws Exception {
		AffymetrixGuidType guid = new AffymetrixGuidType();
		guid.generateGuid();
		int sz = guid.getGuid().length;
		assertEquals(sz, AffymetrixGuidType.GUID_LENGTH);
	}

}
