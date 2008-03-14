/*
 * ArrayDataTest.java
 * JUnit based test
 *
 * Created on October 23, 2005, 7:04 PM
 */

package affymetrix.calvin.array;

import java.util.ArrayList;
import java.util.List;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.array.CreateStep.CreateStepType;
import affymetrix.calvin.parameter.ParameterNameValueDefaultRequired;
import affymetrix.calvin.utils.AffymetrixGuidType;

/**
 * 
 * @author ljevon
 */
public class ArrayDataTest extends TestCase {

	public ArrayDataTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(ArrayDataTest.class);

		return suite;
	}

	/**
	 * Test of getArraySetFileIdentifier method, of class affymetrix.calvin.array.ArrayData.
	 */
	public void testArraySetFileIdentifier() throws Exception {
		ArrayData d = new ArrayData();
		AffymetrixGuidType g = new AffymetrixGuidType();
		g.setGuid("guid");
		d.setArraySetFileIdentifier(g);
		AffymetrixGuidType g2 = d.getArraySetFileIdentifier();
		assertEquals(g2, g);
	}

	/**
	 * Test of getDataTypeIdentifier method, of class affymetrix.calvin.array.ArrayData.
	 */
	public void testDataTypeIdentifier() throws Exception {
		ArrayData d = new ArrayData();
		AffymetrixGuidType g = new AffymetrixGuidType();
		g.setGuid("guid");
		d.setDataTypeIdentifier(g);
		AffymetrixGuidType g2 = d.getDataTypeIdentifier();
		assertEquals(g2, g);
	}

	/**
	 * Test of getCreatedStep method, of class affymetrix.calvin.array.ArrayData.
	 */
	public void testCreatedStep() {
		ArrayData d = new ArrayData();
		d.setCreatedStep(CreateStepType.values()[2]);
		assertEquals(d.getCreatedStep(), 2);
	}

	/**
	 * Test of getInitialProject method, of class affymetrix.calvin.array.ArrayData.
	 */
	public void testInitialProject() {
		ArrayData d = new ArrayData();
		d.setInitialProject("proj");
		assertEquals(d.getInitialProject(), "proj");
	}

	/**
	 * Test of getCreationDateTime method, of class affymetrix.calvin.array.ArrayData.
	 */
	public void testCreationDateTime() {
		ArrayData d = new ArrayData();
		d.setCreationDateTime("date");
		assertEquals(d.getCreationDateTime(), "date");
	}

	/**
	 * Test of getCreatedBy method, of class affymetrix.calvin.array.ArrayData.
	 */
	public void testCreatedBy() {
		ArrayData d = new ArrayData();
		d.setCreatedBy("by");
		assertEquals(d.getCreatedBy(), "by");
	}

	/**
	 * Test of getPhysicalArraysAttributes method, of class affymetrix.calvin.array.ArrayData.
	 */
	public void testPhysicalArraysAttributes() {
		ArrayData d = new ArrayData();
		List<ArrayAttributes> at = new ArrayList<ArrayAttributes>();
		ArrayAttributes attr = new ArrayAttributes();
		at.add(attr);
		d.addPhysicalArraysAttributes(at);
		assertEquals(d.getPhysicalArraysAttributes().get(0), at.get(0));
	}

	/**
	 * Test of getUserAttributes method, of class affymetrix.calvin.array.ArrayData.
	 */
	public void testUserAttributes() {
		ArrayData d = new ArrayData();
		List<ParameterNameValueDefaultRequired> at = new ArrayList<ParameterNameValueDefaultRequired>();
		ParameterNameValueDefaultRequired pnv = new ParameterNameValueDefaultRequired();
		at.add(pnv);
		d.addUserAttributes(at);
		assertEquals(d.getUserAttributes(), at);
	}

}
