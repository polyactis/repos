/*
 * ArrayMediaTest.java
 * JUnit based test
 *
 * Created on October 23, 2005, 7:03 PM
 */

package affymetrix.calvin.array;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.array.ArrayMedia.ArrayMediaType;

/**
 * 
 * @author ljevon
 */
public class ArrayMediaTest extends TestCase {

	public ArrayMediaTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(ArrayMediaTest.class);

		return suite;
	}

	/**
	 * Test of getMedia method, of class affymetrix.calvin.array.ArrayMedia.
	 */
	public void testArrayMedia() {
		assertEquals(ArrayMedia.toString(ArrayMediaType.CartridgeMedia), "Cartridge");
		assertTrue(ArrayMedia.toArrayMediaType("Cartridge") == ArrayMediaType.CartridgeMedia);
		assertTrue(ArrayMedia.toArrayMediaType("PlateOrStrip") == ArrayMediaType.PlateOrStripMedia);
		assertEquals(ArrayMedia.toString(ArrayMediaType.PlateOrStripMedia), "PlateOrStrip");
	}

}
