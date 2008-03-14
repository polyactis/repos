/*
 * CreateStepTest.java
 * JUnit based test
 *
 * Created on October 23, 2005, 7:03 PM
 */

package affymetrix.calvin.array;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.array.CreateStep.CreateStepType;

/**
 * 
 * @author ljevon
 */
public class CreateStepTest extends TestCase {

	public CreateStepTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(CreateStepTest.class);

		return suite;
	}

	/**
	 * Test of getStep method, of class affymetrix.calvin.array.CreateStep.
	 */
	public void testCreateStep() {
		assertTrue(CreateStep.toCreateStepType("ArrayRegistration") == CreateStepType.ArrayRegistrationStep);
		assertEquals(CreateStep.toString(CreateStepType.ArrayRegistrationStep), "ArrayRegistration");
		
		assertTrue(CreateStep.toCreateStepType("CELAnalysis") == CreateStepType.CELAnalysisStep);
		assertEquals(CreateStep.toString(CreateStepType.CELAnalysisStep), "CELAnalysis");
		
		assertTrue(CreateStep.toCreateStepType("Gridding") == CreateStepType.GriddingStep);
		assertEquals(CreateStep.toString(CreateStepType.GriddingStep), "Gridding");
		
		assertTrue(CreateStep.toCreateStepType("None") == CreateStepType.NoStep);
		assertEquals(CreateStep.toString(CreateStepType.NoStep), "None");
		
		assertTrue(CreateStep.toCreateStepType("Other") == CreateStepType.OtherStep);
		assertEquals(CreateStep.toString(CreateStepType.OtherStep), "Other");
	}

}
