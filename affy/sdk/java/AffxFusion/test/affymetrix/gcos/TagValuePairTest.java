/*
 * TagValuePairTest.java
 * JUnit based test
 *
 * Created on October 12, 2005, 4:40 PM
 */

package affymetrix.gcos;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class TagValuePairTest extends TestCase {

	public TagValuePairTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(TagValuePairTest.class);

		return suite;
	}

	/**
	 * Test of getTag method, of class affymetrix.gcos.TagValuePair.
	 */
	public void testTag() {
		TagValuePair t = new TagValuePair();
		t.setTag("tag");
		assertEquals(t.getTag(), "tag");
	}

	/**
	 * Test of getValue method, of class affymetrix.gcos.TagValuePair.
	 */
	public void testValue() {
		TagValuePair t = new TagValuePair();
		t.setValue("value");
		assertEquals(t.getValue(), "value");
	}

	/**
	 * Test of copy method, of class affymetrix.gcos.TagValuePair.
	 */
	public void testCopy() {
		TagValuePair t1 = new TagValuePair();
		t1.setValue("value");
		t1.setTag("tag");
		TagValuePair t2 = new TagValuePair(t1);
		t1.setTag("na");
		t1.setValue("na");
		assertEquals(t2.getTag(), "tag");
		assertEquals(t2.getValue(), "value");

		t1 = t2.copy();
		t2.setTag("na");
		t2.setValue("na");
		assertEquals(t1.getTag(), "tag");
		assertEquals(t1.getValue(), "value");
	}

}
