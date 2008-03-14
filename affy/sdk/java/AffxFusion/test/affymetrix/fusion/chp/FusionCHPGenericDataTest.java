/*
 * FusionCHPGenericDataTest.java
 * JUnit based test
 *
 * Created on December 2, 2005, 12:47 PM
 */

package affymetrix.fusion.chp;

import java.io.File;
import java.util.Arrays;

import junit.framework.TestCase;
import affymetrix.calvin.utils.AffymetrixGuidType;

/**
 * 
 * @author ljevon
 */
public class FusionCHPGenericDataTest extends TestCase {

	public FusionCHPGenericDataTest(String testName) {
		super(testName);
		FusionCHPGenericData.registerReader();
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public void testmethod_FileId() throws Exception {
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_tiling_file");
		FusionCHPData chp = FusionCHPDataReg.read(f.getCanonicalPath());
		assertTrue(chp != null);
		FusionCHPGenericData tileChp = FusionCHPGenericData.fromBase(chp);
		assertTrue(tileChp != null);
		AffymetrixGuidType guid = new AffymetrixGuidType("0000039321-1131034775-0000032391-0000005436-0000004827");
		assertTrue(Arrays.equals(chp.getFileId().getGuid(), guid.getGuid()));
	}

	public void testRead() throws Exception {
		FusionCHPData chp;
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\no_file_exists");
		chp = FusionCHPDataReg.read(f.getCanonicalPath());
		FusionCHPGenericData genchp = FusionCHPGenericData.fromBase(chp);
		assertTrue(chp == null);
		assertTrue(genchp == null);

		f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_genotype_file");
		chp = FusionCHPDataReg.read(f.getCanonicalPath());
		genchp = FusionCHPGenericData.fromBase(chp);
		assertTrue(chp != null);
		assertTrue(genchp != null);

		f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_expression_file");
		chp = FusionCHPDataReg.read(f.getCanonicalPath());
		genchp = FusionCHPGenericData.fromBase(chp);
		assertTrue(chp != null);
		assertTrue(genchp != null);

		f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_tiling_file");
		chp = FusionCHPDataReg.read(f.getCanonicalPath());
		genchp = FusionCHPGenericData.fromBase(chp);
		assertTrue(chp != null);
		assertTrue(genchp != null);

		f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_reseq_file");
		chp = FusionCHPDataReg.read(f.getCanonicalPath());
		genchp = FusionCHPGenericData.fromBase(chp);
		assertTrue(chp != null);
		assertTrue(genchp != null);

		f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_universal_file");
		chp = FusionCHPDataReg.read(f.getCanonicalPath());
		genchp = FusionCHPGenericData.fromBase(chp);
		assertTrue(chp != null);
		assertTrue(genchp != null);

		f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\test.exp.abs.CHP");
		chp = FusionCHPDataReg.read(f.getCanonicalPath());
		genchp = FusionCHPGenericData.fromBase(chp);
		assertTrue(chp == null);
		assertTrue(genchp == null);

		f = new File("..\\..\\..\\..\\calvin_files\\fusion\\data\\small_cel_file_full_datheader");
		chp = FusionCHPDataReg.read(f.getCanonicalPath());
		genchp = FusionCHPGenericData.fromBase(chp);
		assertTrue(chp != null);
		assertTrue(genchp != null);
		AffymetrixGuidType guid = new AffymetrixGuidType("0000065535-1131033797-0000017421-0000012382-0000000292");
		assertTrue(Arrays.equals(chp.getFileId().getGuid(), guid.getGuid()));
		assertTrue(genchp.getData() != null);
	}
}
