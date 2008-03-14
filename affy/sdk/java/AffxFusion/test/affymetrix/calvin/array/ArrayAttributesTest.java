/*
 * ArrayAttributesTest.java
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
import affymetrix.calvin.array.ArrayMedia.ArrayMediaType;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.utils.AffymetrixGuidType;

/**
 * 
 * @author ljevon
 */
public class ArrayAttributesTest extends TestCase {

	public ArrayAttributesTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(ArrayAttributesTest.class);

		return suite;
	}

	/**
	 * Test of getIdentifier method, of class affymetrix.calvin.array.ArrayAttributes.
	 */
	public void testIdentifier() throws Exception {
		ArrayAttributes a = new ArrayAttributes();
		a.setArrayBarcode("bar");
		assertEquals(a.getArrayBarcode(), "bar");
		a.setArrayName("name");
		assertEquals(a.getArrayName(), "name");
		a.setComment("c");
		assertEquals(a.getComment(), "c");
		a.setCreatedBy("by");
		assertEquals(a.getCreatedBy(), "by");
		a.setCreatedStep(CreateStep.CreateStepType.values()[1]);
		assertEquals(a.getCreatedStep(), CreateStep.CreateStepType.values()[1]);
		a.setCreationDateTime("now");
		assertEquals(a.getCreationDateTime(), "now");
		a.setCustomerBarcode("bar2");
		assertEquals(a.getCustomerBarcode(), "bar2");
		AffymetrixGuidType g = new AffymetrixGuidType();
		g.setGuid("guid");
		a.setIdentifier(g);
		AffymetrixGuidType g2 = a.getIdentifier();
		assertTrue(g.equals(g2));
		a.setMasterFile("file");
		assertEquals(a.getMasterFile(), "file");
		a.setMedia(ArrayMediaType.values()[1]);
		assertEquals(a.getMedia(), ArrayMediaType.values()[1]);
		a.setMediaRow(3);
		assertEquals(a.getMediaRow(), 3);
		a.setMediaCol(4);
		assertEquals(a.getMediaCol(), 4);
		a.setPATAssignment(PATAssignmentMethod.PATAssignmentMethodType.values()[3]);
		assertEquals(a.getPATAssignment(), PATAssignmentMethod.PATAssignmentMethodType.values()[3]);
		List<ParameterNameValue> at = new ArrayList<ParameterNameValue>();
		ParameterNameValue pnv = new ParameterNameValue();
		at.add(0, pnv);
		a.getAttributes().addAll(at);
		List<ParameterNameValue> at2 = new ArrayList<ParameterNameValue>();
		at2 = a.getAttributes();
		assertEquals(at, at2);
	}
}
