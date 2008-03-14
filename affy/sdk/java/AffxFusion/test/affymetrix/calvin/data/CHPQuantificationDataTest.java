/*
 * CHPQuantificationDataTest.java
 * JUnit based test
 *
 * Created on December 27, 2005, 8:09 AM
 */

package affymetrix.calvin.data;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class CHPQuantificationDataTest extends TestCase {

	public CHPQuantificationDataTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(CHPQuantificationDataTest.class);

		return suite;
	}

	public void test_FileName() {
		CHPQuantificationData data = new CHPQuantificationData();
		data.setFilename("file");
		assertEquals(data.getFilename(), "file");
	}

	public void test_ArrayType() {
		CHPQuantificationData data = new CHPQuantificationData();
		assertEquals(data.getArrayType(), null);
	}

	public void test_AlgName() {
		CHPQuantificationData data = new CHPQuantificationData();
		assertEquals(data.getAlgName(), null);
	}

	public void test_AlgVersion() {
		CHPQuantificationData data = new CHPQuantificationData();
		assertEquals(data.getAlgVersion(), null);
	}

}
