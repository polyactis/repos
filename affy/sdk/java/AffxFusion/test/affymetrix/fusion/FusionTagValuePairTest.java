/*
 * FusionTagValuePairTest.java
 * JUnit based test
 *
 * Created on December 5, 2005, 9:35 AM
 */

package affymetrix.fusion;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.parameter.ParameterNameValue;

/**
 * 
 * @author ljevon
 */
public class FusionTagValuePairTest extends TestCase {

	public FusionTagValuePairTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(FusionTagValuePairTest.class);

		return suite;
	}

	/**
	 * Test of getTag method, of class affymetrix.fusion.FusionTagValuePair.
	 */
	public void testTag() {
		FusionTagValuePair p = new FusionTagValuePair();
		p.setTag("tag");
		assertEquals(p.getTag(), "tag");
	}

	/**
	 * Test of getValue method, of class affymetrix.fusion.FusionTagValuePair.
	 */
	public void testValue() {
		FusionTagValuePair p = new FusionTagValuePair();
		p.setValue("value");
		assertEquals(p.getValue(), "value");
	}

	/**
	 * Test of getDetailed method, of class affymetrix.fusion.FusionTagValuePair.
	 */
	public void testDetailed() {
		FusionTagValuePair p = new FusionTagValuePair();
		ParameterNameValue det = new ParameterNameValue();
		det.setName("tag");
		det.setValueText("value");
		p.setDetailed(det);
		det = null;
		det = p.getDetailed();
		assertEquals(det.getName(), "tag");
		assertEquals(det.getValueText(), "value");
	}

	/**
	 * Test of copy method, of class affymetrix.fusion.FusionTagValuePair.
	 */
	public void testCopy() {
		FusionTagValuePair p = new FusionTagValuePair();
		p.setTag("tag");
		p.setValue("value");
		FusionTagValuePair pc = p.copy();
		assertEquals(pc.getTag(), "tag");
		assertEquals(pc.getValue(), "value");
		assertEquals(pc.getDetailed(), null);

		ParameterNameValue det = new ParameterNameValue();
		det.setName("tag");
		det.setValueText("value");
		p.setDetailed(det);
		pc = p.copy();
		assertEquals(pc.getTag(), "tag");
		assertEquals(pc.getValue(), "value");
		assertEquals(pc.getDetailed().getName(), "tag");
		assertEquals(pc.getDetailed().getValueText(), "value");

	}

}
