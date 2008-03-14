/*
 * PATAssignmentMethodTest.java
 * JUnit based test
 *
 * Created on October 23, 2005, 7:03 PM
 */

package affymetrix.calvin.array;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.array.PATAssignmentMethod.PATAssignmentMethodType;

/**
 * 
 * @author ljevon
 */
public class PATAssignmentMethodTest extends TestCase {

	public PATAssignmentMethodTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(PATAssignmentMethodTest.class);

		return suite;
	}

	/**
	 * Test of getMethod method, of class affymetrix.calvin.array.PATAssignmentMethod.
	 */
	public void testMethod() {
		assertTrue(PATAssignmentMethod.toPATAssignmentMethodType(PATAssignmentMethod.PAT_ASSIGNMENT_BARCODE) == PATAssignmentMethodType.AffyBarcodeAssignment);
		assertTrue(PATAssignmentMethod.toPATAssignmentMethodType(PATAssignmentMethod.PAT_ASSIGNMENT_NONE) == PATAssignmentMethodType.NoAssignment);
		assertTrue(PATAssignmentMethod.toPATAssignmentMethodType(PATAssignmentMethod.PAT_ASSIGNMENT_OTHER) == PATAssignmentMethodType.OtherAssignment);
		assertTrue(PATAssignmentMethod.toPATAssignmentMethodType(PATAssignmentMethod.PAT_ASSIGNMENT_USER_SELECTED) == PATAssignmentMethodType.UserSelectedAssignment);
	}

	/**
	 * Test of toString method, of class affymetrix.calvin.array.PATAssignmentMethod.
	 */
	public void testToString() {
		assertEquals(PATAssignmentMethod.toString(PATAssignmentMethodType.AffyBarcodeAssignment), PATAssignmentMethod.PAT_ASSIGNMENT_BARCODE);
		assertEquals(PATAssignmentMethod.toString(PATAssignmentMethodType.NoAssignment), PATAssignmentMethod.PAT_ASSIGNMENT_NONE);
		assertEquals(PATAssignmentMethod.toString(PATAssignmentMethodType.OtherAssignment), PATAssignmentMethod.PAT_ASSIGNMENT_OTHER);
		assertEquals(PATAssignmentMethod.toString(PATAssignmentMethodType.UserSelectedAssignment), PATAssignmentMethod.PAT_ASSIGNMENT_USER_SELECTED);
	}

}
