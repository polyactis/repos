/*
 * EXPFileDataTest.java
 * JUnit based test
 *
 * Created on October 20, 2005, 8:38 AM
 */

package affymetrix.gcos.exp;

import java.io.File;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.gcos.TagValuePair;

/**
 * 
 * @author ljevon
 */
public class EXPFileDataTest extends TestCase {

	public EXPFileDataTest(String testName) throws Exception {
		super(testName);
		FULL_EXP = new File("..\\..\\..\\CPPTest\\data\\full.EXP").getCanonicalPath();
		BARE_EXP = new File("..\\..\\..\\CPPTest\\data\\bare.EXP").getCanonicalPath();
		ABORT_EXP = new File("..\\..\\..\\CPPTest\\data\\hyb_abort.EXP").getCanonicalPath();
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(EXPFileDataTest.class);

		return suite;
	}

	private String FULL_EXP;

	private String BARE_EXP;

	private String ABORT_EXP;

	public void testproperty_FileName() {
		EXPFileData exp = new EXPFileData();
		String path = "test";
		exp.setFileName(path);
		assertEquals(path, exp.getFileName());
	}

	public void testproperty_ArrayType() {
		EXPFileData exp = new EXPFileData();
		String path = FULL_EXP;
		exp.setFileName(path);
		assertEquals(exp.exists(), true);
		exp.read();
		assertEquals(exp.getArrayType(), "Test3");

		exp.clear();
		assertEquals(exp.getArrayType(), "");
	}

	public void testmethod_exists() {
		EXPFileData exp = new EXPFileData();
		String path = FULL_EXP;
		exp.setFileName(path);
		assertEquals(exp.exists(), true);
	}

	public void testmethod_existsWhenFileNotexists() {
		EXPFileData exp = new EXPFileData();
		String path = "no_file_exists";
		exp.setFileName(path);
		assertEquals(exp.exists(), false);
	}

	public void testmethod_read_full() {
		EXPFileData exp = new EXPFileData();
		String path = FULL_EXP;
		exp.setFileName(path);
		assertEquals(exp.exists(), true);
		assertEquals(exp.read(), true);

		assertEquals(exp.getArrayType(), "Test3");

		TagValuePair param;
		assertEquals(exp.getNumScanParameters(), 7);
		param = exp.getScanParameter(0);
		assertEquals(param.getTag(), "Pixel Size");
		assertEquals(param.getValue(), "3");

		param = exp.getScanParameter(1);
		assertEquals(param.getTag(), "Filter");
		assertEquals(param.getValue(), "570");

		param = exp.getScanParameter(2);
		assertEquals(param.getTag(), "Scan Temperature");
		assertEquals(param.getValue(), "200");

		param = exp.getScanParameter(3);
		assertEquals(param.getTag(), "Scan Date");
		assertEquals(param.getValue(), "Aug 01 2002 08:37AM");

		param = exp.getScanParameter(4);
		assertEquals(param.getTag(), "Scanner ID");
		assertEquals(param.getValue(), "123123");

		param = exp.getScanParameter(5);
		assertEquals(param.getTag(), "Number of Scans");
		assertEquals(param.getValue(), "1");

		param = exp.getScanParameter(6);
		assertEquals(param.getTag(), "Scanner Type");
		assertEquals(param.getValue(), "HP");

		assertEquals(exp.getNumHybParameters(), 19);
		param = exp.getHybParameter(0);
		assertEquals(param.getTag(), "Protocol");
		assertEquals(param.getValue(), "EukGE-WS1v4");

		param = exp.getHybParameter(1);
		assertEquals(param.getTag(), "Wash A1 Recovery Mixes");
		assertEquals(param.getValue(), "0");

		param = exp.getHybParameter(2);
		assertEquals(param.getTag(), "Wash A1 Temperature (C)");
		assertEquals(param.getValue(), "25");

		param = exp.getHybParameter(3);
		assertEquals(param.getTag(), "Number of Wash A1 Cycles");
		assertEquals(param.getValue(), "10");

		param = exp.getHybParameter(4);
		assertEquals(param.getTag(), "Mixes per Wash A1 Cycle");
		assertEquals(param.getValue(), "2");

		param = exp.getHybParameter(5);
		assertEquals(param.getTag(), "Wash B Recovery Mixes");
		assertEquals(param.getValue(), "0");

		param = exp.getHybParameter(6);
		assertEquals(param.getTag(), "Wash B Temperature (C)");
		assertEquals(param.getValue(), "50");

		param = exp.getHybParameter(7);
		assertEquals(param.getTag(), "Number of Wash B Cycles");
		assertEquals(param.getValue(), "4");

		param = exp.getHybParameter(8);
		assertEquals(param.getTag(), "Mixes per Wash B Cycle");
		assertEquals(param.getValue(), "15");

		param = exp.getHybParameter(9);
		assertEquals(param.getTag(), "Stain Temperature (C)");
		assertEquals(param.getValue(), "25");

		param = exp.getHybParameter(10);
		assertEquals(param.getTag(), "Stain Time (seconds)");
		assertEquals(param.getValue(), "1800");

		param = exp.getHybParameter(11);
		assertEquals(param.getTag(), "Wash A2 Recovery Mixes");
		assertEquals(param.getValue(), "0");

		param = exp.getHybParameter(12);
		assertEquals(param.getTag(), "Wash A2 Temperature (C)");
		assertEquals(param.getValue(), "25");

		param = exp.getHybParameter(13);
		assertEquals(param.getTag(), "Number of Wash A2 Cycles");
		assertEquals(param.getValue(), "10");

		param = exp.getHybParameter(14);
		assertEquals(param.getTag(), "Mixes per Wash A2 Cycle");
		assertEquals(param.getValue(), "4");

		param = exp.getHybParameter(15);
		assertEquals(param.getTag(), "Holding Temperature (C)");
		assertEquals(param.getValue(), "25");

		param = exp.getHybParameter(16);
		assertEquals(param.getTag(), "Station");
		assertEquals(param.getValue(), "1");

		param = exp.getHybParameter(17);
		assertEquals(param.getTag(), "Module");
		assertEquals(param.getValue(), "1");

		param = exp.getHybParameter(18);
		assertEquals(param.getTag(), "Hybridize Date");
		assertEquals(param.getValue(), "Oct 08 2002 03:31PM");

		assertEquals(exp.getNumSampleParameters(), 8);
		param = exp.getSampleParameter(0);
		assertEquals(param.getTag(), "Chip Lot");
		assertEquals(param.getValue(), "123");

		param = exp.getSampleParameter(1);
		assertEquals(param.getTag(), "Operator");
		assertEquals(param.getValue(), "ljevon");

		param = exp.getSampleParameter(2);
		assertEquals(param.getTag(), "Sample Type");
		assertEquals(param.getValue(), "s_type");

		param = exp.getSampleParameter(3);
		assertEquals(param.getTag(), "Description");
		assertEquals(param.getValue(), "Demo data");

		param = exp.getSampleParameter(4);
		assertEquals(param.getTag(), "Project");
		assertEquals(param.getValue(), "s_proj");

		param = exp.getSampleParameter(5);
		assertEquals(param.getTag(), "Comments");
		assertEquals(param.getValue(), "None");

		param = exp.getSampleParameter(6);
		assertEquals(param.getTag(), "Solution Type");
		assertEquals(param.getValue(), "sol type");

		param = exp.getSampleParameter(7);
		assertEquals(param.getTag(), "Solution Lot");
		assertEquals(param.getValue(), "123123");

	}

	public void testmethod_read_bare() {
		EXPFileData exp = new EXPFileData();
		String path = BARE_EXP;
		exp.setFileName(path);
		assertEquals(exp.exists(), true);
		assertEquals(exp.read(), true);

		assertEquals(exp.getArrayType(), "Mapping10K_Xba131");

		TagValuePair param;
		assertEquals(exp.getNumScanParameters(), 1);
		param = exp.getScanParameter(0);
		assertEquals(param.getTag(), "Scanner Type");
		assertEquals(param.getValue(), "M18");

		assertEquals(exp.getNumHybParameters(), 0);

		assertEquals(exp.getNumSampleParameters(), 1);
		param = exp.getSampleParameter(0);
		assertEquals(param.getTag(), "Operator");
		assertEquals(param.getValue(), "sgelle");
	}

	public void testmethod_read_hyb_abort() {
		EXPFileData exp = new EXPFileData();
		String path = ABORT_EXP;
		exp.setFileName(path);
		assertEquals(exp.exists(), true);
		assertEquals(exp.read(), true);

		assertEquals(exp.getArrayType(), "HG-U133A");

		TagValuePair param;
		assertEquals(exp.getNumScanParameters(), 0);

		assertEquals(exp.getNumHybParameters(), 7);

		param = exp.getHybParameter(0);
		assertEquals(param.getTag(), "Protocol");
		assertEquals(param.getValue(), "EukGE-WS1v4");

		param = exp.getHybParameter(1);
		assertEquals(param.getTag(), "REMOVE VIAL");
		assertEquals(param.getValue(), "");

		param = exp.getHybParameter(2);
		assertEquals(param.getTag(), "            20°C");
		assertEquals(param.getValue(), "");

		param = exp.getHybParameter(3);
		assertEquals(param.getTag(), "Aborted by user");
		assertEquals(param.getValue(), "");

		param = exp.getHybParameter(4);
		assertEquals(param.getTag(), "Station");
		assertEquals(param.getValue(), "1");

		param = exp.getHybParameter(5);
		assertEquals(param.getTag(), "Module");
		assertEquals(param.getValue(), "2");

		param = exp.getHybParameter(6);
		assertEquals(param.getTag(), "Hybridize Date");
		assertEquals(param.getValue(), "Mar 28 2003 09:52AM");

		assertEquals(exp.getNumSampleParameters(), 0);
	}

	public void testmethod_readWhenFileNotexists() {
		EXPFileData exp = new EXPFileData();
		String path = "test";
		exp.setFileName(path);
		assertEquals(exp.read(), false);
	}

}
