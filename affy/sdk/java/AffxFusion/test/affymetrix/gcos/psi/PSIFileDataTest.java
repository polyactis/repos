/*
 * PSIFileDataTest.java
 * JUnit based test
 *
 * Created on October 12, 2005, 3:45 PM
 */

package affymetrix.gcos.psi;

import java.io.File;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class PSIFileDataTest extends TestCase {

	public PSIFileDataTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(PSIFileDataTest.class);

		return suite;
	}

	/**
	 * Test of getFileName method, of class affymetrix.gcos.psi.PSIFileData.
	 */
	public void testFileName() {
		PSIFileData ps = new PSIFileData();
		ps.setFileName("name");
		assertEquals(ps.getFileName(), "name");
	}

	/**
	 * Test of GetProbeSetCount method, of class affymetrix.gcos.psi.PSIFileData.
	 */
	public void testGetProbeSetCount() {
		PSIFileData ps = new PSIFileData();
		assertEquals(ps.getProbeSetCount(), 0);
	}

	/**
	 * Test of Read method, of class affymetrix.gcos.psi.PSIFileData.
	 */
	public void testRead() throws Exception {
		PSIFileData ps = new PSIFileData();
		ps.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.PSI").getCanonicalPath());
		ps.read();
		assertEquals(ps.getProbeSetCount(), 5);
		assertEquals(ps.getProbeSetName(0), "one");
		assertEquals(ps.getProbeSetName(1), "two");
		assertEquals(ps.getProbeSetName(2), "three");
		assertEquals(ps.getProbeSetName(3), "four");
		assertEquals(ps.getProbeSetName(4), "five");
		assertEquals(ps.getProbePairs(0), 1);
		assertEquals(ps.getProbePairs(1), 2);
		assertEquals(ps.getProbePairs(2), 3);
		assertEquals(ps.getProbePairs(3), 4);
		assertEquals(ps.getProbePairs(4), 5);

	}

	/**
	 * Test of Exists method, of class affymetrix.gcos.psi.PSIFileData.
	 */
	public void testExists() throws Exception {
		PSIFileData ps = new PSIFileData();
		ps.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.PSI").getCanonicalPath());
		assertEquals(ps.exists(), true);

		ps.setFileName("C:\\cvs_head\\affy\\sdk\\file\\CPPTest\\data\\nofile.PSI");
		assertEquals(ps.exists(), false);

	}

	/**
	 * Test of Clear method, of class affymetrix.gcos.psi.PSIFileData.
	 */
	public void testClear() throws Exception {
		PSIFileData ps = new PSIFileData();
		ps.setFileName(new File("..\\..\\..\\CPPTest\\data\\test.PSI").getCanonicalPath());
		ps.read();
		assertEquals(ps.getProbeSetCount(), 5);
		ps.clear();
		assertEquals(ps.getProbeSetCount(), 0);
	}

}
