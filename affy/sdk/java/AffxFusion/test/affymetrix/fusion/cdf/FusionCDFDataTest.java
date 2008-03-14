/*
 * FusionCDFDataTest.java
 * JUnit based test
 *
 * Created on October 28, 2005, 7:45 AM
 */

package affymetrix.fusion.cdf;

import java.io.File;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class FusionCDFDataTest extends TestCase {

	public FusionCDFDataTest(String testName) throws Exception {
		super(testName);
		// use this if data is mounted
		EXP_XDA_FILE = new File("\\\\affymetrix.com\\shares\\santaclara\\TestDataFiles\\CDF\\Test3-xda.CDF")
				.getCanonicalPath();
		EXP_ASCII_FILE = new File("\\\\affymetrix.com\\shares\\santaclara\\TestDataFiles\\CDF\\Test3-ascii.CDF")
				.getCanonicalPath();
		NO_FILE = new File("\\\\affymetrix.com\\shares\\santaclara\\TestDataFiles\\CDF\\NoFile.CDF").getCanonicalPath();
		MAP_XDA_FILE = new File("\\\\affymetrix.com\\shares\\santaclara\\TestDataFiles\\CDF\\Mapping10K_Xba131-xda.CDF")
				.getCanonicalPath();
		MAP_ASCII_FILE = new File("\\\\affymetrix.com\\shares\\santaclara\\TestDataFiles\\CDF\\Mapping10K_Xba131-ascii.CDF")
				.getCanonicalPath();

		// use this if data is local
		EXP_XDA_FILE = new File("D:\\TestDataFiles\\CDF\\Test3-xda.CDF").getCanonicalPath();
		EXP_ASCII_FILE = new File("D:\\TestDataFiles\\CDF\\Test3-ascii.CDF").getCanonicalPath();
		NO_FILE = new File("D:\\TestDataFiles\\CDF\\NoFile.CDF").getCanonicalPath();
		MAP_XDA_FILE = new File("D:\\TestDataFiles\\CDF\\Mapping10K_Xba131-xda.CDF").getCanonicalPath();
		MAP_ASCII_FILE = new File("D:\\TestDataFiles\\CDF\\Mapping10K_Xba131-ascii.CDF").getCanonicalPath();
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(FusionCDFDataTest.class);

		return suite;
	}

	private String EXP_XDA_FILE;

	private String EXP_ASCII_FILE;

	private String NO_FILE;

	private String MAP_XDA_FILE;

	private String MAP_ASCII_FILE;

	public void test_FusionCDFHeader() {
		FusionCDFHeader header = new FusionCDFHeader();
		assertEquals(header.getCols(), 0);
		assertEquals(header.getRows(), 0);
		assertEquals(header.getNumProbeSets(), 0);
		assertEquals(header.getNumQCProbeSets(), 0);
		assertEquals(header.getReference(), null);
	}

	public void test_FusionCDFProbeInformation() {
		FusionCDFProbeInformation probe = new FusionCDFProbeInformation();
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 0);
		assertEquals(probe.getY(), 0);
		assertEquals(probe.getPBase(), ' ');
		assertEquals(probe.getTBase(), ' ');
	}

	public void test_FusionCDFProbeGroupInformation() {
		FusionCDFProbeGroupInformation group = new FusionCDFProbeGroupInformation();
		assertEquals(group.getDirection(), FusionDirectionType.NoDirection);
		assertEquals(group.getNumLists(), 0);
		assertEquals(group.getNumCells(), 0);
		assertEquals(group.getNumCellsPerList(), 0);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), 0);
		assertEquals(group.getName(), null);
	}

	public void test_FusionCDFProbeSetInformation() {
		FusionCDFProbeSetInformation set = new FusionCDFProbeSetInformation();
		assertEquals(set.getProbeSetType(), FusionGeneChipProbeSetType.UnknownProbeSetType);
		assertEquals(set.getDirection(), FusionDirectionType.NoDirection);
		assertEquals(set.getNumLists(), 0);
		assertEquals(set.getNumGroups(), 0);
		assertEquals(set.getNumCells(), 0);
		assertEquals(set.getNumCellsPerList(), 0);
		assertEquals(set.getProbeSetNumber(), 0);
	}

	public void test_FusionCDFQCProbeInformation() {
		FusionCDFQCProbeInformation probe = new FusionCDFQCProbeInformation();
		assertEquals(probe.getX(), 0);
		assertEquals(probe.getY(), 0);
		assertEquals(probe.getProbeLength(), 0);
		assertFalse(probe.isPMProbe());
		assertFalse(probe.isBackground());
	}

	public void test_CDFQCProbeSetInformation() {
		FusionCDFQCProbeSetInformation set = new FusionCDFQCProbeSetInformation();
		assertEquals(set.getNumCells(), 0);
		assertEquals(set.getQCProbeSetType(), FusionGeneChipQCProbeSetType.UnknownQCProbeSetType);
	}

	public void testmethod_IsXDACompatibleFile() {
		FusionCDFData cdf = new FusionCDFData();
		cdf.setFileName(EXP_XDA_FILE);
		assertTrue(cdf.isXDACompatibleFile());
		cdf.setFileName(EXP_ASCII_FILE);
		assertFalse(cdf.isXDACompatibleFile());
	}

	public void testmethod_Exists() {
		FusionCDFData cdf = new FusionCDFData();
		cdf.setFileName(EXP_XDA_FILE);
		assertTrue(cdf.exists());
		cdf.setFileName(NO_FILE);
		assertFalse(cdf.exists());
	}

	public void testmethod_ReadHeader_with_ASCII() {
		FusionCDFData cdf = new FusionCDFData();
		cdf.setFileName(EXP_ASCII_FILE);
		assertTrue(cdf.readHeader());
		assertEquals(cdf.getChipType(), "Test3-ascii");
		FusionCDFHeader header = cdf.getHeader();
		assertEquals(header.getCols(), 126);
		assertEquals(header.getRows(), 126);
		assertEquals(header.getNumProbeSets(), 345);
		assertEquals(header.getNumQCProbeSets(), 13);
		assertEquals(header.getReference(), null);
	}

	public void testmethod_ReadHeader_with_XDA() {
		FusionCDFData cdf = new FusionCDFData();
		cdf.setFileName(EXP_XDA_FILE);
		assertTrue(cdf.readHeader());
		assertEquals(cdf.getChipType(), "Test3-xda");
		FusionCDFHeader header = cdf.getHeader();
		assertEquals(header.getCols(), 126);
		assertEquals(header.getRows(), 126);
		assertEquals(header.getNumProbeSets(), 345);
		assertEquals(header.getNumQCProbeSets(), 13);
		assertEquals(header.getReference(), null);
	}

	public void test_ExpressionXDA() {
		FusionCDFData cdf = new FusionCDFData();
		cdf.setFileName(EXP_XDA_FILE);
		assertTrue(cdf.read());

		assertEquals(cdf.getChipType(), "Test3-xda");
		FusionCDFHeader header = cdf.getHeader();
		assertEquals(header.getCols(), 126);
		assertEquals(header.getRows(), 126);
		assertEquals(header.getNumProbeSets(), 345);
		assertEquals(header.getNumQCProbeSets(), 13);
		assertEquals(header.getReference(), null);

		assertEquals(cdf.getProbeSetName(0), "Pae_16SrRNA_s_at");
		assertEquals(cdf.getProbeSetName(1), "Pae_23SrRNA_s_at");
		assertEquals(cdf.getProbeSetName(2), "PA1178_oprH_at");
		assertEquals(cdf.getProbeSetName(header.getNumProbeSets() - 1), "AFFX_ratb2/X14115_at");

		for (int i = 0; i < header.getNumProbeSets(); i++) {
			assertEquals(cdf.getProbeSetType(i), FusionGeneChipProbeSetType.ExpressionProbeSetType);
		}

		FusionCDFProbeSetInformation set = new FusionCDFProbeSetInformation();
		FusionCDFProbeGroupInformation group = new FusionCDFProbeGroupInformation();
		FusionCDFProbeInformation probe = new FusionCDFProbeInformation();

		cdf.getProbeSetInformation(0, set);
		assertEquals(set.getProbeSetType(), FusionGeneChipProbeSetType.ExpressionProbeSetType);
		assertEquals(set.getDirection(), FusionDirectionType.AntiSenseDirection);
		assertEquals(set.getNumLists(), 16);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 32);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 1000);

		set.getGroup(0, group);
		assertEquals(group.getDirection(), FusionDirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 16);
		assertEquals(group.getNumCells(), 32);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "Pae_16SrRNA_s_at");

		group.getCell(0, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 111);
		assertEquals(probe.getY(), 79);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		group.getCell(1, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 111);
		assertEquals(probe.getY(), 80);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		group.getCell(group.getNumCells() - 1, probe);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), group.getNumLists() - 1);
		assertEquals(probe.getX(), 93);
		assertEquals(probe.getY(), 82);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("g"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("g"), 0);

		cdf.getProbeSetInformation(1, set);
		assertEquals(set.getProbeSetType(), FusionGeneChipProbeSetType.ExpressionProbeSetType);
		assertEquals(set.getDirection(), FusionDirectionType.AntiSenseDirection);
		assertEquals(set.getNumLists(), 16);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 32);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 1001);

		set.getGroup(0, group);
		assertEquals(group.getDirection(), FusionDirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 16);
		assertEquals(group.getNumCells(), 32);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "Pae_23SrRNA_s_at");

		group.getCell(0, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 124);
		assertEquals(probe.getY(), 95);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("t"), 0);
		group.getCell(1, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 124);
		assertEquals(probe.getY(), 96);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("t"), 0);
		group.getCell(group.getNumCells() - 1, probe);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), group.getNumLists() - 1);
		assertEquals(probe.getX(), 36);
		assertEquals(probe.getY(), 8);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		cdf.getProbeSetInformation(header.getNumProbeSets() - 1, set);
		assertEquals(set.getProbeSetType(), FusionGeneChipProbeSetType.ExpressionProbeSetType);
		assertEquals(set.getDirection(), FusionDirectionType.AntiSenseDirection);
		assertEquals(set.getNumLists(), 20);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 40);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 3101);

		set.getGroup(0, group);
		assertEquals(group.getDirection(), FusionDirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 20);
		assertEquals(group.getNumCells(), 40);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "AFFX_ratb2/X14115_at");

		group.getCell(0, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 20);
		assertEquals(probe.getY(), 113);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		group.getCell(1, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 20);
		assertEquals(probe.getY(), 114);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		group.getCell(group.getNumCells() - 1, probe);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), group.getNumLists() - 1);
		assertEquals(probe.getX(), 119);
		assertEquals(probe.getY(), 58);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		FusionCDFQCProbeSetInformation qcset = new FusionCDFQCProbeSetInformation();
		FusionCDFQCProbeInformation qcprobe = new FusionCDFQCProbeInformation();

		cdf.getQCProbeSetInformation(0, qcset);
		assertEquals(qcset.getNumCells(), 300);
		assertEquals(qcset.getQCProbeSetType(), FusionGeneChipQCProbeSetType.GeneExpNegativeQCProbeSetType);

		qcset.getCell(0, qcprobe);
		assertEquals(qcprobe.getX(), 77);
		assertEquals(qcprobe.getY(), 82);
		assertEquals(qcprobe.getProbeLength(), 20);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcset.getCell(1, qcprobe);
		assertEquals(qcprobe.getX(), 77);
		assertEquals(qcprobe.getY(), 83);
		assertEquals(qcprobe.getProbeLength(), 20);
		assertTrue(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcset.getCell(qcset.getNumCells() - 1, qcprobe);
		assertEquals(qcprobe.getX(), 15);
		assertEquals(qcprobe.getY(), 86);
		assertEquals(qcprobe.getProbeLength(), 1);
		assertFalse(qcprobe.isPMProbe());
		assertTrue(qcprobe.isBackground());

		cdf.getQCProbeSetInformation(header.getNumQCProbeSets() - 1, qcset);
		assertEquals(qcset.getNumCells(), 9);
		assertEquals(qcset.getQCProbeSetType(), FusionGeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);

		qcset.getCell(0, qcprobe);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 60);
		assertEquals(qcprobe.getProbeLength(), 25);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcset.getCell(1, qcprobe);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 61);
		assertEquals(qcprobe.getProbeLength(), 25);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcset.getCell(qcset.getNumCells() - 1, qcprobe);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 64);
		assertEquals(qcprobe.getProbeLength(), 25);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		cdf.getQCProbeSetInformationByType(FusionGeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType, qcset);
		assertEquals(qcset.getNumCells(), 9);
		assertEquals(qcset.getQCProbeSetType(), FusionGeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);

		qcset.getCell(0, qcprobe);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 60);
		assertEquals(qcprobe.getProbeLength(), 25);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcset.getCell(1, qcprobe);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 61);
		assertEquals(qcprobe.getProbeLength(), 25);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcset.getCell(qcset.getNumCells() - 1, qcprobe);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 64);
		assertEquals(qcprobe.getProbeLength(), 25);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());
	}

	public void test_ExpressionASCII() {
		FusionCDFData cdf = new FusionCDFData();
		cdf.setFileName(EXP_ASCII_FILE);
		assertTrue(cdf.read());

		assertEquals(cdf.getChipType(), "Test3-ascii");
		FusionCDFHeader header = cdf.getHeader();
		assertEquals(header.getCols(), 126);
		assertEquals(header.getRows(), 126);
		assertEquals(header.getNumProbeSets(), 345);
		assertEquals(header.getNumQCProbeSets(), 13);
		assertEquals(header.getReference(), null);

		assertEquals(cdf.getProbeSetName(0), "Pae_16SrRNA_s_at");
		assertEquals(cdf.getProbeSetName(1), "Pae_23SrRNA_s_at");
		assertEquals(cdf.getProbeSetName(2), "PA1178_oprH_at");
		assertEquals(cdf.getProbeSetName(header.getNumProbeSets() - 1), "AFFX_ratb2/X14115_at");

		for (int i = 0; i < header.getNumProbeSets(); i++) {
			assertEquals(cdf.getProbeSetType(i), FusionGeneChipProbeSetType.ExpressionProbeSetType);
		}

		FusionCDFProbeSetInformation set = new FusionCDFProbeSetInformation();
		FusionCDFProbeGroupInformation group = new FusionCDFProbeGroupInformation();
		FusionCDFProbeInformation probe = new FusionCDFProbeInformation();

		cdf.getProbeSetInformation(0, set);
		assertEquals(set.getProbeSetType(), FusionGeneChipProbeSetType.ExpressionProbeSetType);
		assertEquals(set.getDirection(), FusionDirectionType.AntiSenseDirection);
		assertEquals(set.getNumLists(), 16);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 32);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 1000);

		set.getGroup(0, group);
		assertEquals(group.getDirection(), FusionDirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 16);
		assertEquals(group.getNumCells(), 32);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "Pae_16SrRNA_s_at");

		group.getCell(0, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 111);
		assertEquals(probe.getY(), 79);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		group.getCell(1, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 111);
		assertEquals(probe.getY(), 80);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		group.getCell(group.getNumCells() - 1, probe);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), group.getNumLists() - 1);
		assertEquals(probe.getX(), 93);
		assertEquals(probe.getY(), 82);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("g"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("g"), 0);

		cdf.getProbeSetInformation(1, set);
		assertEquals(set.getProbeSetType(), FusionGeneChipProbeSetType.ExpressionProbeSetType);
		assertEquals(set.getDirection(), FusionDirectionType.AntiSenseDirection);
		assertEquals(set.getNumLists(), 16);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 32);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 1001);

		set.getGroup(0, group);
		assertEquals(group.getDirection(), FusionDirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 16);
		assertEquals(group.getNumCells(), 32);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "Pae_23SrRNA_s_at");

		group.getCell(0, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 124);
		assertEquals(probe.getY(), 95);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("t"), 0);
		group.getCell(1, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 124);
		assertEquals(probe.getY(), 96);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("t"), 0);
		group.getCell(group.getNumCells() - 1, probe);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), group.getNumLists() - 1);
		assertEquals(probe.getX(), 36);
		assertEquals(probe.getY(), 8);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		cdf.getProbeSetInformation(header.getNumProbeSets() - 1, set);
		assertEquals(set.getProbeSetType(), FusionGeneChipProbeSetType.ExpressionProbeSetType);
		assertEquals(set.getDirection(), FusionDirectionType.AntiSenseDirection);
		assertEquals(set.getNumLists(), 20);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 40);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 3101);

		set.getGroup(0, group);
		assertEquals(group.getDirection(), FusionDirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 20);
		assertEquals(group.getNumCells(), 40);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "AFFX_ratb2/X14115_at");

		group.getCell(0, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 20);
		assertEquals(probe.getY(), 113);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		group.getCell(1, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 20);
		assertEquals(probe.getY(), 114);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		group.getCell(group.getNumCells() - 1, probe);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), group.getNumLists() - 1);
		assertEquals(probe.getX(), 119);
		assertEquals(probe.getY(), 58);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		FusionCDFQCProbeSetInformation qcset = new FusionCDFQCProbeSetInformation();
		FusionCDFQCProbeInformation qcprobe = new FusionCDFQCProbeInformation();

		cdf.getQCProbeSetInformation(0, qcset);
		assertEquals(qcset.getNumCells(), 300);
		assertEquals(qcset.getQCProbeSetType(), FusionGeneChipQCProbeSetType.GeneExpNegativeQCProbeSetType);

		qcset.getCell(0, qcprobe);
		assertEquals(qcprobe.getX(), 77);
		assertEquals(qcprobe.getY(), 82);
		assertEquals(qcprobe.getProbeLength(), 20);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcset.getCell(1, qcprobe);
		assertEquals(qcprobe.getX(), 77);
		assertEquals(qcprobe.getY(), 83);
		assertEquals(qcprobe.getProbeLength(), 20);

		qcset.getCell(qcset.getNumCells() - 1, qcprobe);
		assertEquals(qcprobe.getX(), 15);
		assertEquals(qcprobe.getY(), 86);
		assertEquals(qcprobe.getProbeLength(), 1);

		cdf.getQCProbeSetInformation(header.getNumQCProbeSets() - 1, qcset);
		assertEquals(qcset.getNumCells(), 9);
		assertEquals(qcset.getQCProbeSetType(), FusionGeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);

		qcset.getCell(0, qcprobe);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 60);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcset.getCell(1, qcprobe);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 61);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcset.getCell(qcset.getNumCells() - 1, qcprobe);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 64);
		assertEquals(qcprobe.getProbeLength(), 25);

		cdf.getQCProbeSetInformationByType(FusionGeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType, qcset);
		assertEquals(qcset.getNumCells(), 9);
		assertEquals(qcset.getQCProbeSetType(), FusionGeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);

		qcset.getCell(0, qcprobe);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 60);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcset.getCell(1, qcprobe);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 61);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcset.getCell(qcset.getNumCells() - 1, qcprobe);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 64);
		assertEquals(qcprobe.getProbeLength(), 25);
	}

	public void test_GenotypingXDA() {
		FusionCDFData cdf = new FusionCDFData();
		cdf.setFileName(MAP_XDA_FILE);
		assertTrue(cdf.read());

		assertEquals(cdf.getChipType(), "Mapping10K_Xba131-xda");
		FusionCDFHeader header = cdf.getHeader();
		assertEquals(header.getCols(), 712);
		assertEquals(header.getRows(), 712);
		assertEquals(header.getNumProbeSets(), 11564);
		assertEquals(header.getNumQCProbeSets(), 9);
		assertEquals(header.getReference(), null);

		assertEquals(cdf.getProbeSetName(0), "AFFX-5Q-123");
		assertEquals(cdf.getProbeSetName(1), "AFFX-5Q-456");
		assertEquals(cdf.getProbeSetName(2), "AFFX-5Q-789");
		assertEquals(cdf.getProbeSetName(header.getNumProbeSets() - 1), "SNP_A-1508078");

		for (int i = 0; i < header.getNumProbeSets(); i++) {
			assertEquals(cdf.getProbeSetType(i), FusionGeneChipProbeSetType.GenotypingProbeSetType);
		}

		FusionCDFProbeSetInformation set = new FusionCDFProbeSetInformation();
		FusionCDFProbeGroupInformation group = new FusionCDFProbeGroupInformation();
		FusionCDFProbeInformation probe = new FusionCDFProbeInformation();

		cdf.getProbeSetInformation(0, set);
		assertEquals(set.getProbeSetType(), FusionGeneChipProbeSetType.GenotypingProbeSetType);
		assertEquals(set.getDirection(), FusionDirectionType.SenseDirection);
		assertEquals(set.getNumLists(), 30);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 60);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 41);

		set.getGroup(0, group);
		assertEquals(group.getDirection(), FusionDirectionType.SenseDirection);
		assertEquals(group.getNumLists(), 30);
		assertEquals(group.getNumCells(), 60);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "AFFX-5Q-123");

		group.getCell(0, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 12);
		assertEquals(probe.getX(), 323);
		assertEquals(probe.getY(), 386);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		group.getCell(1, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 12);
		assertEquals(probe.getX(), 323);
		assertEquals(probe.getY(), 387);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		group.getCell(group.getNumCells() - 1, probe);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), 46);
		assertEquals(probe.getX(), 337);
		assertEquals(probe.getY(), 389);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		cdf.getProbeSetInformation(header.getNumProbeSets() - 1, set);
		assertEquals(set.getProbeSetType(), FusionGeneChipProbeSetType.GenotypingProbeSetType);
		assertEquals(set.getDirection(), FusionDirectionType.SenseDirection);
		assertEquals(set.getNumLists(), 20);
		assertEquals(set.getNumGroups(), 4);
		assertEquals(set.getNumCells(), 40);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 12548);

		set.getGroup(0, group);
		assertEquals(group.getDirection(), FusionDirectionType.SenseDirection);
		assertEquals(group.getNumLists(), 5);
		assertEquals(group.getNumCells(), 10);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), 4);
		assertEquals(group.getName(), "SNP_A-1508078A");

		group.getCell(0, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 642);
		assertEquals(probe.getY(), 665);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("g"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("c"), 0);
		group.getCell(1, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 642);
		assertEquals(probe.getY(), 666);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("c"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("c"), 0);
		group.getCell(group.getNumCells() - 1, probe);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), 21);
		assertEquals(probe.getX(), 172);
		assertEquals(probe.getY(), 699);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("t"), 0);

		set.getGroup(3, group);
		assertEquals(group.getDirection(), FusionDirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 5);
		assertEquals(group.getNumCells(), 10);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 15);
		assertEquals(group.getStop(), 19);
		assertEquals(group.getName(), "SNP_A-1508078G");

		group.getCell(0, probe);
		assertEquals(probe.getListIndex(), 15);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 648);
		assertEquals(probe.getY(), 604);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("g"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("g"), 0);
		group.getCell(1, probe);
		assertEquals(probe.getListIndex(), 15);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 648);
		assertEquals(probe.getY(), 603);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("c"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("g"), 0);
		group.getCell(group.getNumCells() - 1, probe);
		assertEquals(probe.getListIndex(), 19);
		assertEquals(probe.getExpos(), 21);
		assertEquals(probe.getX(), 220);
		assertEquals(probe.getY(), 686);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		FusionCDFQCProbeSetInformation qcset = new FusionCDFQCProbeSetInformation();
		FusionCDFQCProbeInformation qcprobe = new FusionCDFQCProbeInformation();

		cdf.getQCProbeSetInformation(0, qcset);
		assertEquals(qcset.getNumCells(), 24);
		assertEquals(qcset.getQCProbeSetType(), FusionGeneChipQCProbeSetType.HybPositiveQCProbeSetType);

		qcset.getCell(0, qcprobe);
		assertEquals(qcprobe.getX(), 473);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcset.getCell(1, qcprobe);
		assertEquals(qcprobe.getX(), 236);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcset.getCell(qcset.getNumCells() - 1, qcprobe);
		assertEquals(qcprobe.getX(), 711);
		assertEquals(qcprobe.getY(), 472);
		assertEquals(qcprobe.getProbeLength(), 16);

		cdf.getQCProbeSetInformation(header.getNumQCProbeSets() - 1, qcset);
		assertEquals(qcset.getNumCells(), 9);
		assertEquals(qcset.getQCProbeSetType(), FusionGeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);

		qcset.getCell(0, qcprobe);
		assertEquals(qcprobe.getX(), 356);
		assertEquals(qcprobe.getY(), 354);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcset.getCell(1, qcprobe);
		assertEquals(qcprobe.getX(), 356);
		assertEquals(qcprobe.getY(), 355);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcset.getCell(qcset.getNumCells() - 1, qcprobe);
		assertEquals(qcprobe.getX(), 356);
		assertEquals(qcprobe.getY(), 358);
		assertEquals(qcprobe.getProbeLength(), 25);

		cdf.getQCProbeSetInformationByType(FusionGeneChipQCProbeSetType.HybPositiveQCProbeSetType, qcset);
		assertEquals(qcset.getNumCells(), 24);
		assertEquals(qcset.getQCProbeSetType(), FusionGeneChipQCProbeSetType.HybPositiveQCProbeSetType);

		qcset.getCell(0, qcprobe);
		assertEquals(qcprobe.getX(), 473);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcset.getCell(1, qcprobe);
		assertEquals(qcprobe.getX(), 236);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcset.getCell(qcset.getNumCells() - 1, qcprobe);
		assertEquals(qcprobe.getX(), 711);
		assertEquals(qcprobe.getY(), 472);
		assertEquals(qcprobe.getProbeLength(), 16);

	}

	public void test_GenotypingASCII() {
		FusionCDFData cdf = new FusionCDFData();
		cdf.setFileName(MAP_ASCII_FILE);
		assertTrue(cdf.read());

		assertEquals(cdf.getChipType(), "Mapping10K_Xba131-ascii");
		FusionCDFHeader header = cdf.getHeader();
		assertEquals(header.getCols(), 712);
		assertEquals(header.getRows(), 712);
		assertEquals(header.getNumProbeSets(), 11564);
		assertEquals(header.getNumQCProbeSets(), 9);
		assertEquals(header.getReference(), null);

		assertEquals(cdf.getProbeSetName(0), "AFFX-5Q-123");
		assertEquals(cdf.getProbeSetName(1), "AFFX-5Q-456");
		assertEquals(cdf.getProbeSetName(2), "AFFX-5Q-789");
		assertEquals(cdf.getProbeSetName(header.getNumProbeSets() - 1), "SNP_A-1508078");

		for (int i = 0; i < header.getNumProbeSets(); i++) {
			assertEquals(cdf.getProbeSetType(i), FusionGeneChipProbeSetType.GenotypingProbeSetType);
		}

		FusionCDFProbeSetInformation set = new FusionCDFProbeSetInformation();
		FusionCDFProbeGroupInformation group = new FusionCDFProbeGroupInformation();
		FusionCDFProbeInformation probe = new FusionCDFProbeInformation();

		cdf.getProbeSetInformation(0, set);
		assertEquals(set.getProbeSetType(), FusionGeneChipProbeSetType.GenotypingProbeSetType);
		assertEquals(set.getDirection(), FusionDirectionType.SenseDirection);
		assertEquals(set.getNumLists(), 30);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 60);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 41);

		set.getGroup(0, group);
		assertEquals(group.getDirection(), FusionDirectionType.SenseDirection);
		assertEquals(group.getNumLists(), 30);
		assertEquals(group.getNumCells(), 60);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "AFFX-5Q-123");

		group.getCell(0, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 12);
		assertEquals(probe.getX(), 323);
		assertEquals(probe.getY(), 386);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		group.getCell(1, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 12);
		assertEquals(probe.getX(), 323);
		assertEquals(probe.getY(), 387);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		group.getCell(group.getNumCells() - 1, probe);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), 46);
		assertEquals(probe.getX(), 337);
		assertEquals(probe.getY(), 389);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		cdf.getProbeSetInformation(header.getNumProbeSets() - 1, set);
		assertEquals(set.getProbeSetType(), FusionGeneChipProbeSetType.GenotypingProbeSetType);
		assertEquals(set.getDirection(), FusionDirectionType.SenseDirection);
		assertEquals(set.getNumLists(), 20);
		assertEquals(set.getNumGroups(), 4);
		assertEquals(set.getNumCells(), 40);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 12548);

		set.getGroup(0, group);
		assertEquals(group.getDirection(), FusionDirectionType.SenseDirection);
		assertEquals(group.getNumLists(), 5);
		assertEquals(group.getNumCells(), 10);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), 4);
		assertEquals(group.getName(), "SNP_A-1508078A");

		group.getCell(0, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 642);
		assertEquals(probe.getY(), 665);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("g"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("c"), 0);
		group.getCell(1, probe);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 642);
		assertEquals(probe.getY(), 666);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("c"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("c"), 0);
		group.getCell(group.getNumCells() - 1, probe);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), 21);
		assertEquals(probe.getX(), 172);
		assertEquals(probe.getY(), 699);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("t"), 0);

		set.getGroup(3, group);
		assertEquals(group.getDirection(), FusionDirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 5);
		assertEquals(group.getNumCells(), 10);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 15);
		assertEquals(group.getStop(), 19);
		assertEquals(group.getName(), "SNP_A-1508078G");

		group.getCell(0, probe);
		assertEquals(probe.getListIndex(), 15);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 648);
		assertEquals(probe.getY(), 604);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("g"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("g"), 0);
		group.getCell(1, probe);
		assertEquals(probe.getListIndex(), 15);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 648);
		assertEquals(probe.getY(), 603);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("c"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("g"), 0);
		group.getCell(group.getNumCells() - 1, probe);
		assertEquals(probe.getListIndex(), 19);
		assertEquals(probe.getExpos(), 21);
		assertEquals(probe.getX(), 220);
		assertEquals(probe.getY(), 686);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		FusionCDFQCProbeSetInformation qcset = new FusionCDFQCProbeSetInformation();
		FusionCDFQCProbeInformation qcprobe = new FusionCDFQCProbeInformation();

		cdf.getQCProbeSetInformation(0, qcset);
		assertEquals(qcset.getNumCells(), 24);
		assertEquals(qcset.getQCProbeSetType(), FusionGeneChipQCProbeSetType.HybPositiveQCProbeSetType);

		qcset.getCell(0, qcprobe);
		assertEquals(qcprobe.getX(), 473);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcset.getCell(1, qcprobe);
		assertEquals(qcprobe.getX(), 236);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcset.getCell(qcset.getNumCells() - 1, qcprobe);
		assertEquals(qcprobe.getX(), 711);
		assertEquals(qcprobe.getY(), 472);
		assertEquals(qcprobe.getProbeLength(), 16);

		cdf.getQCProbeSetInformation(header.getNumQCProbeSets() - 1, qcset);
		assertEquals(qcset.getNumCells(), 9);
		assertEquals(qcset.getQCProbeSetType(), FusionGeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);

		qcset.getCell(0, qcprobe);
		assertEquals(qcprobe.getX(), 356);
		assertEquals(qcprobe.getY(), 354);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcset.getCell(1, qcprobe);
		assertEquals(qcprobe.getX(), 356);
		assertEquals(qcprobe.getY(), 355);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcset.getCell(qcset.getNumCells() - 1, qcprobe);
		assertEquals(qcprobe.getX(), 356);
		assertEquals(qcprobe.getY(), 358);
		assertEquals(qcprobe.getProbeLength(), 25);

		cdf.getQCProbeSetInformationByType(FusionGeneChipQCProbeSetType.HybPositiveQCProbeSetType, qcset);
		assertEquals(qcset.getNumCells(), 24);
		assertEquals(qcset.getQCProbeSetType(), FusionGeneChipQCProbeSetType.HybPositiveQCProbeSetType);

		qcset.getCell(0, qcprobe);
		assertEquals(qcprobe.getX(), 473);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcset.getCell(1, qcprobe);
		assertEquals(qcprobe.getX(), 236);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcset.getCell(qcset.getNumCells() - 1, qcprobe);
		assertEquals(qcprobe.getX(), 711);
		assertEquals(qcprobe.getY(), 472);
		assertEquals(qcprobe.getProbeLength(), 16);
	}

}
