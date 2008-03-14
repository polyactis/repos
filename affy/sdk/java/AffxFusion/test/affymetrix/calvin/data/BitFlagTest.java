/*
 * BitFlagTest.java
 * JUnit based test
 *
 * Created on October 23, 2005, 7:05 PM
 */

package affymetrix.calvin.data;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class BitFlagTest extends TestCase {

	public BitFlagTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(BitFlagTest.class);

		return suite;
	}

	/**
	 * Test of clear method, of class affymetrix.calvin.data.BitFlag.
	 */
	public void testClass() {
		BitFlag flags = new BitFlag();
		assertEquals(flags.getFlags(), 0);
		flags.setFlags((short)27);
		assertEquals(flags.getFlags(), 27);
		flags.clear();
		assertEquals(flags.getFlags(), 0);
		flags.setDefaultDataSetHdr(true);
		assertEquals(flags.hasDefaultDataSetHdr(), true);
		BitFlag flags2 = new BitFlag((short)33);
		assertEquals(flags2.getFlags(), 33);
		BitFlag flags3 = new BitFlag();
		flags3 = flags2;
		assertEquals(flags3.getFlags(), 33);
	}

}
