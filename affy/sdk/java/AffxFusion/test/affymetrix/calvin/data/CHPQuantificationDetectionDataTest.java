/*
 * CHPQuantificationDetectionDataTest.java
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
public class CHPQuantificationDetectionDataTest extends TestCase {

	public CHPQuantificationDetectionDataTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(CHPQuantificationDetectionDataTest.class);

		return suite;
	}

	public void test_FileName() {
		CHPQuantificationDetectionData data = new CHPQuantificationDetectionData();
		data.setFilename("file");
		assertEquals(data.getFilename(), "file");
	}

	public void test_ArrayType() {
		CHPQuantificationDetectionData data = new CHPQuantificationDetectionData();
		assertEquals(data.getArrayType(), null);
	}

	public void test_AlgName() {
		CHPQuantificationDetectionData data = new CHPQuantificationDetectionData();
		assertEquals(data.getAlgName(), null);
	}

	public void test_AlgVersion() {
		CHPQuantificationDetectionData data = new CHPQuantificationDetectionData();
		assertEquals(data.getAlgVersion(), null);
	}

}
