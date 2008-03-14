/*
 * GenericDataHeaderTest.java
 * JUnit based test
 *
 * Created on October 23, 2005, 7:05 PM
 */

package affymetrix.calvin.data;

import java.io.UnsupportedEncodingException;
import java.util.Arrays;
import java.util.List;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.utils.AffymetrixGuidType;

/**
 * 
 * @author ljevon
 */
public class GenericDataHeaderTest extends TestCase {

	public GenericDataHeaderTest(String testName) {
		super(testName);
	}

	private GenericDataHeader header;

	protected void setUp() throws Exception {
		header = new GenericDataHeader();
	}

	protected void tearDown() throws Exception {
		header = null;
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(GenericDataHeaderTest.class);

		return suite;
	}

	public void testFileTypeIdTest() {
		String p1 = "file type ID";
		header.setFileTypeId(p1);
		String p2 = header.getFileTypeId();
		assertEquals(p1, p2);
	}

	public void testFileIdTest() throws Exception {
		AffymetrixGuidType p1 = new AffymetrixGuidType();
		p1.generateGuid();
		header.setFileId(p1);
		byte[] p2 = header.getFileId().getGuid();
		assertTrue(Arrays.equals(p1.getGuid(), p2));
	}

	public void testFileCreationTimeTest() {
		String p1 = "20050301T12:42:11Z";
		header.setFileCreationTime(p1);
		String p2 = header.getFileCreationTime();
		assertEquals(p1, p2);
	}

	public void testLocaleTest() {
		String p1 = "locale";
		header.setLocale(p1);
		String p2 = header.getLocale();
		assertEquals(p1, p2);
	}

	public void testNameValTest() {
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

		List<ParameterNameValue> params = header.getNameValParams();
		assertEquals(params.size(), 3);

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

	public void testAddParentEntryTest() throws UnsupportedEncodingException {
		GenericDataHeader header1 = new GenericDataHeader();
		AffymetrixGuidType g = new AffymetrixGuidType();
		byte[] s1 = g.generateGuid();
		header1.setFileId(g);
		GenericDataHeader header2 = new GenericDataHeader();
		g = new AffymetrixGuidType();
		byte[] s2 = g.generateGuid();
		header2.setFileId(g);
		GenericDataHeader header3 = new GenericDataHeader();
		g = new AffymetrixGuidType();
		byte[] s3 = g.generateGuid();
		header3.setFileId(g);
		header.addParent(header1);
		header.addParent(header2);
		header.addParent(header3);

		List<GenericDataHeader> parents = header.getParents();
		GenericDataHeader p = parents.get(0);
		assertTrue(Arrays.equals(p.getFileId().getGuid(), s1));
		p = header.getParent(0);
		assertTrue(Arrays.equals(p.getFileId().getGuid(), s1));
		p = parents.get(1);
		assertTrue(Arrays.equals(p.getFileId().getGuid(), s2));
		p = header.getParent(1);
		assertTrue(Arrays.equals(p.getFileId().getGuid(), s2));
		p = parents.get(2);
		assertTrue(Arrays.equals(p.getFileId().getGuid(), s3));
		p = header.getParent(2);
		assertTrue(Arrays.equals(p.getFileId().getGuid(), s3));
	}

	public void testGetNumParentsTest() {
		GenericDataHeader header1 = new GenericDataHeader();
		GenericDataHeader header2 = new GenericDataHeader();
		GenericDataHeader header3 = new GenericDataHeader();
		header.addParent(header1);
		header.addParent(header2);
		header.addParent(header3);
		int sz = header.getParentCnt();
		assertEquals(sz, 3);
	}

	public void testGetNameValPairCntTest() {
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
	}

	public void testGetNameValPairTest() {
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
		assertEquals(header.getNameValParam(0).getName(), p[0]);
		assertEquals(header.getNameValParam(0).getValueText(), p[1]);
		assertEquals(header.getNameValParam(1).getName(), p[2]);
		assertEquals(header.getNameValParam(1).getValueText(), p[3]);
		assertEquals(header.getNameValParam(2).getName(), p[4]);
		assertEquals(header.getNameValParam(2).getValueText(), p[5]);
	}

	public void testUpdateNameValPairTest() {
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
		assertEquals(header.getNameValParam(0).getValueText(), p[1]);
	}

	public void testFindNameValPairTest() {
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

		ParameterNameValue t4 = header.findNameValParam(p[0]);
		assertEquals(t4.getName(), p[0]);
		assertEquals(t4.getValueText(), p[1]);
		t4 = header.findNameValParam("Oingo Boingo");
		assertEquals(t4, null);
	}

	public void testGetNameValParamsBeginsWithTest() {
		String names[] = { "prefix-one", "two", "prefix-three" };
		ParameterNameValue t;

		for (int i = 0; i < 3; ++i) {
			t = new ParameterNameValue();
			t.setName(names[i]);
			t.setValueInt32(i);
			header.addNameValParam(t);
		}

		List<ParameterNameValue> v = header.getNameValParamsBeginsWith("prefix-");
		assertEquals(v.size(), 2);

		t = v.get(0);
		assertEquals(t.getName(), names[0]);

		t = v.get(1);
		assertEquals(t.getName(), names[2]);
	}

	public void testFindParentByFileTypeIdTest() {
		GenericDataHeader header1 = new GenericDataHeader();
		header1.setFileTypeId("Fred");
		GenericDataHeader header2 = new GenericDataHeader();
		header2.setFileTypeId("Wilma");
		GenericDataHeader header3 = new GenericDataHeader();
		header3.setFileTypeId("Betty");
		header.addParent(header1);
		header.addParent(header2);
		header.addParent(header3);

		GenericDataHeader gdh = header.findParent("Wilma");
		assertTrue(gdh != null);
		assertEquals(gdh.getFileTypeId(), "Wilma");
		gdh = header.findParent("nonexistant");
		assertEquals(gdh, null);
	}

}
