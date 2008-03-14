/*
 * DataSetHeaderTest.java
 * JUnit based test
 *
 * Created on October 23, 2005, 7:05 PM
 */

package affymetrix.calvin.data;

import java.util.List;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.data.ColumnInfo.DataSetColumnTypes;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.portability.UInt;

/**
 * 
 * @author ljevon
 */
public class DataSetHeaderTest extends TestCase {

	public DataSetHeaderTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(DataSetHeaderTest.class);

		return suite;
	}

	public void testNameTest() {
		DataSetHeader header = new DataSetHeader();
		String p1 = "some_name";
		header.setName(p1);
		String p2 = header.getName();
		assertTrue(p1 == p2);
	}

	public void testNameValPairCntTest() {
		DataSetHeader header = new DataSetHeader();
		String p[] = { "Uma", "Thurman", "Angelina", "Jolie", "Isabella ", "Rossellini" };
		ParameterNameValue t1 = new ParameterNameValue();
		t1.setName(p[0]);
		t1.setValueText(p[1]);
		header.addNameValParam(t1);

		ParameterNameValue t2 = new ParameterNameValue();
		t2.setName(p[2]);
		t2.setValueText(p[3]);
		header.addNameValParam(t2);

		ParameterNameValue t3 = new ParameterNameValue();
		t3.setName(p[4]);
		t3.setValueText(p[5]);
		header.addNameValParam(t3);
		assertTrue(header.getNameValParamCnt() == 3);
	}

	public void testNameValPairTest() {
		DataSetHeader header = new DataSetHeader();
		String p[] = { "Uma", "Thurman", "Angelina", "Jolie", "Isabella ", "Rossellini" };
		ParameterNameValue t1 = new ParameterNameValue();
		t1.setName(p[0]);
		t1.setValueText(p[1]);
		header.addNameValParam(t1);

		ParameterNameValue t2 = new ParameterNameValue();
		t2.setName(p[2]);
		t2.setValueText(p[3]);
		header.addNameValParam(t2);

		ParameterNameValue t3 = new ParameterNameValue();
		t3.setName(p[4]);
		t3.setValueText(p[5]);
		header.addNameValParam(t3);

		assertEquals(header.getNameValParamCnt(), 3);
		List<ParameterNameValue> params = header.getNameValParameters();

		ParameterNameValue param = params.get(0);
		assertEquals(param.getName(), p[0]);
		assertEquals(param.getValueText(), p[1]);

		param = params.get(1);
		assertEquals(param.getName(), p[2]);
		assertEquals(param.getValueText(), p[3]);

		param = params.get(2);
		assertEquals(param.getName(), p[4]);
		assertEquals(param.getValueText(), p[5]);
	}

	public void testColumnNameTest() {
		DataSetHeader header = new DataSetHeader();
		header.addByteColumn("Fred");
		header.addIntColumn("Wilma");
		header.addShortColumn("");
		assertEquals(header.getColumnInfo(0).getName(), "Fred");
		assertEquals(header.getColumnInfo(1).getName(), "Wilma");
		assertEquals(header.getColumnInfo(2).getName(), "");
	}

	public void testColumnValueTypeTest() {
		DataSetHeader header = new DataSetHeader();
		header.addByteColumn("");
		header.addIntColumn("");
		header.addShortColumn("");
		assertTrue(header.getColumnInfo(0).getColumnType() == DataSetColumnTypes.ByteColType);
		assertTrue(header.getColumnInfo(1).getColumnType() == DataSetColumnTypes.IntColType);
		assertTrue(header.getColumnInfo(2).getColumnType() == DataSetColumnTypes.ShortColType);
	}

	public void testColumnSizeTest() {
		DataSetHeader header = new DataSetHeader();
		header.addByteColumn("");
		header.addIntColumn("");
		header.addShortColumn("");
		assertTrue(header.getColumnInfo(0).getSize() == 1);
		assertTrue(header.getColumnInfo(1).getSize() == 4);
		assertTrue(header.getColumnInfo(2).getSize() == 2);
	}

	public void testRowCntTest() {
		DataSetHeader header = new DataSetHeader();
		int rows = 6;
		header.setRowCnt(rows);
		assertTrue(header.getRowCnt() == rows);
	}

	public void testDataStartFilePosTest() throws UnsignedOutOfLimitsException {
		DataSetHeader header = new DataSetHeader();
		assertTrue(header.getDataStartFilePos().equals(0));
		header.setDataStartFilePos(new UInt(12));
		assertTrue(header.getDataStartFilePos().equals(12));
	}

	public void testHeaderStartFilePosTest() throws UnsignedOutOfLimitsException {
		DataSetHeader header = new DataSetHeader();
		assertTrue(header.getHeaderStartFilePos().equals(0));
		header.setHeaderStartFilePos(new UInt(61209));
		assertTrue(header.getHeaderStartFilePos().equals(61209));
	}

}
