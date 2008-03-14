/*
 * GenotypeProbeSetResultsTest.java
 * JUnit based test
 *
 * Created on October 13, 2005, 9:20 AM
 */

package affymetrix.gcos.chp;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class GenotypeProbeSetResultsTest extends TestCase {

	public GenotypeProbeSetResultsTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(GenotypeProbeSetResultsTest.class);

		return suite;
	}

	/**
	 * Test of getAlleleCall method, of class affymetrix.gcos.chp.GenotypeProbeSetResults.
	 */
	public void testAlleleCall() {
		GenotypeProbeSetResults g = new GenotypeProbeSetResults();
		g.setAlleleCall(GenotypeProbeSetResults.ALLELE_AB_CALL);
		assertEquals(g.getAlleleCall(), GenotypeProbeSetResults.ALLELE_AB_CALL);
	}

	/**
	 * Test of setAlleleCall method, of class affymetrix.gcos.chp.GenotypeProbeSetResults.
	 */
	public void testConfidence() {
		GenotypeProbeSetResults g = new GenotypeProbeSetResults();
		assertEquals(g.getConfidence(), 0.0f, 0.00000001f);
		g.setConfidence(1.0f);
		assertEquals(g.getConfidence(), 1.0f, 0.00000001f);
	}

	/**
	 * Test of getRAS1 method, of class affymetrix.gcos.chp.GenotypeProbeSetResults.
	 */
	public void testRAS1() {
		GenotypeProbeSetResults g = new GenotypeProbeSetResults();
		assertEquals(g.getRAS1(), 0.0f, 0.00000001f);
		g.setRAS1(1.0f);
		assertEquals(g.getRAS1(), 1.0f, 0.00000001f);
	}

	/**
	 * Test of getRAS2 method, of class affymetrix.gcos.chp.GenotypeProbeSetResults.
	 */
	public void testRAS2() {
		GenotypeProbeSetResults g = new GenotypeProbeSetResults();
		assertEquals(g.getRAS2(), 0.0f, 0.00000001f);
		g.setRAS2(1.0f);
		assertEquals(g.getRAS2(), 1.0f, 0.00000001f);
	}

	/**
	 * Test of getPValue_AA method, of class affymetrix.gcos.chp.GenotypeProbeSetResults.
	 */
	public void testPValue_AA() {
		GenotypeProbeSetResults g = new GenotypeProbeSetResults();
		assertEquals(g.getPValue_AA(), 0.0f, 0.00000001f);
		g.setPValue_AA(1.0f);
		assertEquals(g.getPValue_AA(), 1.0f, 0.00000001f);
	}

	/**
	 * Test of getPValue_AB method, of class affymetrix.gcos.chp.GenotypeProbeSetResults.
	 */
	public void testPValue_AB() {
		GenotypeProbeSetResults g = new GenotypeProbeSetResults();
		assertEquals(g.getPValue_AB(), 0.0f, 0.00000001f);
		g.setPValue_AB(1.0f);
		assertEquals(g.getPValue_AB(), 1.0f, 0.00000001f);
	}

	/**
	 * Test of getPValue_BB method, of class affymetrix.gcos.chp.GenotypeProbeSetResults.
	 */
	public void testPValue_BB() {
		GenotypeProbeSetResults g = new GenotypeProbeSetResults();
		assertEquals(g.getPValue_BB(), 0.0f, 0.00000001f);
		g.setPValue_BB(1.0f);
		assertEquals(g.getPValue_BB(), 1.0f, 0.00000001f);
	}

	/**
	 * Test of getPValue_NoCall method, of class affymetrix.gcos.chp.GenotypeProbeSetResults.
	 */
	public void testPValue_NoCall() {
		GenotypeProbeSetResults g = new GenotypeProbeSetResults();
		assertEquals(g.getPValue_NoCall(), 0.0f, 0.00000001f);
		g.setPValue_NoCall(1.0f);
		assertEquals(g.getPValue_NoCall(), 1.0f, 0.00000001f);
	}

	/**
	 * Test of getAlleleCallString method, of class affymetrix.gcos.chp.GenotypeProbeSetResults.
	 */
	public void testGetAlleleCallString() {
		GenotypeProbeSetResults g = new GenotypeProbeSetResults();
		g.setAlleleCall(GenotypeProbeSetResults.ALLELE_A_CALL);
		assertEquals(g.getAlleleCallString(), "A");
		g.setAlleleCall(GenotypeProbeSetResults.ALLELE_B_CALL);
		assertEquals(g.getAlleleCallString(), "B");
		g.setAlleleCall(GenotypeProbeSetResults.ALLELE_AB_CALL);
		assertEquals(g.getAlleleCallString(), "AB");
		g.setAlleleCall(GenotypeProbeSetResults.ALLELE_NO_CALL);
		assertEquals(g.getAlleleCallString(), "No Call");
	}

}
