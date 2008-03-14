/*
 * ColumnInfoTest.java
 * JUnit based test
 *
 * Created on October 23, 2005, 7:05 PM
 */

package affymetrix.calvin.data;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.data.ColumnInfo.DataSetColumnTypes;

/**
 * 
 * @author ljevon
 */
public class ColumnInfoTest extends TestCase {

	public ColumnInfoTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(ColumnInfoTest.class);

		return suite;
	}

	/**
	 * Test of equals method, of class affymetrix.calvin.data.ColumnInfo.
	 */
	public void testBase() {
		ColumnInfo c = new ColumnInfo("name", DataSetColumnTypes.ByteColType, 1, 2, 3);
		assertEquals(c.getName(), "name");
		assertEquals(c.getColumnType(), DataSetColumnTypes.ByteColType);
		assertEquals(c.getLength(), 2);
		assertEquals(c.getSize(), 5);
	}

}
