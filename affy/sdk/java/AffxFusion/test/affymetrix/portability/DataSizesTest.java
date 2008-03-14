/*
 * DataSizesTest.java
 * JUnit based test
 *
 * Created on October 25, 2005, 8:07 PM
 */

package affymetrix.portability;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class DataSizesTest extends TestCase {

	public DataSizesTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(DataSizesTest.class);

		return suite;
	}

	public void testSizes() {
		assertEquals(DataSizes.CHAR_SIZE, 1);
		assertEquals(DataSizes.INT_SIZE, 4);
		assertEquals(DataSizes.FLOAT_SIZE, 4);
		assertEquals(DataSizes.CHAR_SIZE, 1);
	}

}
