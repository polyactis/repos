/*
 * CDFFileDataTest.java
 * JUnit based test
 *
 * Created on October 19, 2005, 8:18 AM
 */

package affymetrix.gcos.cdf;

import java.io.File;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class CDFFileDataTest extends TestCase {

	public CDFFileDataTest(String testName) throws Exception {
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
		TestSuite suite = new TestSuite(CDFFileDataTest.class);

		return suite;
	}

	private String EXP_XDA_FILE;

	private String EXP_ASCII_FILE;

	private String NO_FILE;

	private String MAP_XDA_FILE;

	private String MAP_ASCII_FILE;

	public void test_CDFFileHeader() {
		CDFFileHeader header = new CDFFileHeader();
		assertEquals(header.getCols(), 0);
		assertEquals(header.getRows(), 0);
		assertEquals(header.getNumProbeSets(), 0);
		assertEquals(header.getNumQCProbeSets(), 0);
		assertEquals(header.getReference(), null);
	}

	public void test_CDFProbeInformation() {
		CDFProbeInformation probe = new CDFProbeInformation();
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 0);
		assertEquals(probe.getY(), 0);
		assertEquals(probe.getPBase(), ' ');
		assertEquals(probe.getTBase(), ' ');
	}

	public void test_CDFProbeGroupInformation() {
		CDFProbeGroupInformation group = new CDFProbeGroupInformation();
		assertEquals(group.getDirection(), DirectionType.NoDirection);
		assertEquals(group.getNumLists(), 0);
		assertEquals(group.getNumCells(), 0);
		assertEquals(group.getNumCellsPerList(), 0);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), 0);
		assertEquals(group.getName(), "");
	}

	public void test_CDFProbeSetInformation() {
		CDFProbeSetInformation set = new CDFProbeSetInformation();
		assertEquals(set.getProbeSetType(), GeneChipProbeSetType.UnknownProbeSetType);
		assertEquals(set.getDirection(), DirectionType.NoDirection);
		assertEquals(set.getNumLists(), 0);
		assertEquals(set.getNumGroups(), 0);
		assertEquals(set.getNumCells(), 0);
		assertEquals(set.getNumCellsPerList(), 0);
		assertEquals(set.getProbeSetNumber(), 0);
	}

	public void test_CDFQCProbeInformation() {
		CDFQCProbeInformation probe = new CDFQCProbeInformation();
		assertEquals(probe.getX(), 0);
		assertEquals(probe.getY(), 0);
		assertEquals(probe.getProbeLength(), 0);
		assertFalse(probe.isPMProbe());
		assertFalse(probe.isBackground());
	}

	public void test_CDFQCProbeSetInformation() {
		CDFQCProbeSetInformation set = new CDFQCProbeSetInformation();
		assertEquals(set.getNumCells(), 0);
		assertEquals(set.getQCProbeSetType(), GeneChipQCProbeSetType.UnknownQCProbeSetType);
	}

	public void testmethod_IsXDACompatibleFile() {
		CDFFileData cdf = new CDFFileData();
		cdf.setFileName(EXP_XDA_FILE);
		assertTrue(cdf.isXDACompatibleFile());
		cdf.setFileName(EXP_ASCII_FILE);
		assertFalse(cdf.isXDACompatibleFile());
	}

	public void testmethod_Exists() {
		CDFFileData cdf = new CDFFileData();
		cdf.setFileName(EXP_XDA_FILE);
		assertTrue(cdf.exists());
		cdf.setFileName(NO_FILE);
		assertFalse(cdf.exists());
	}

	public void testmethod_ReadHeader_with_ASCII() {
		CDFFileData cdf = new CDFFileData();
		cdf.setFileName(EXP_ASCII_FILE);
		assertTrue(cdf.readHeader());
		assertEquals(cdf.getChipType(), "Test3-ascii");
		CDFFileHeader header = cdf.getHeader();
		assertEquals(header.getCols(), 126);
		assertEquals(header.getRows(), 126);
		assertEquals(header.getNumProbeSets(), 345);
		assertEquals(header.getNumQCProbeSets(), 13);
		assertEquals(header.getReference(), null);
	}

	public void testmethod_ReadHeader_with_XDA() {
		CDFFileData cdf = new CDFFileData();
		cdf.setFileName(EXP_XDA_FILE);
		assertTrue(cdf.readHeader());
		assertEquals(cdf.getChipType(), "Test3-xda");
		CDFFileHeader header = cdf.getHeader();
		assertEquals(header.getCols(), 126);
		assertEquals(header.getRows(), 126);
		assertEquals(header.getNumProbeSets(), 345);
		assertEquals(header.getNumQCProbeSets(), 13);
		assertEquals(header.getReference(), null);
	}

	public void test_ExpressionXDA() {
		CDFFileData cdf = new CDFFileData();
		cdf.setFileName(EXP_XDA_FILE);
		assertTrue(cdf.read());

		assertEquals(cdf.getChipType(), "Test3-xda");
		CDFFileHeader header = cdf.getHeader();
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
			assertEquals(cdf.getProbeSetType(i), GeneChipProbeSetType.ExpressionProbeSetType);
		}

		CDFProbeSetInformation set;
		CDFProbeGroupInformation group;
		CDFProbeInformation probe;

		set = cdf.getProbeSetInformation(0);
		assertEquals(set.getProbeSetType(), GeneChipProbeSetType.ExpressionProbeSetType);
		assertEquals(set.getDirection(), DirectionType.AntiSenseDirection);
		assertEquals(set.getNumLists(), 16);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 32);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 1000);

		group = set.getGroup(0);
		assertEquals(group.getDirection(), DirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 16);
		assertEquals(group.getNumCells(), 32);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "Pae_16SrRNA_s_at");

		probe = group.getCell(0);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 111);
		assertEquals(probe.getY(), 79);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		probe = group.getCell(1);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 111);
		assertEquals(probe.getY(), 80);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		probe = group.getCell(group.getNumCells() - 1);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), group.getNumLists() - 1);
		assertEquals(probe.getX(), 93);
		assertEquals(probe.getY(), 82);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("g"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("g"), 0);

		set = cdf.getProbeSetInformation(1);
		assertEquals(set.getProbeSetType(), GeneChipProbeSetType.ExpressionProbeSetType);
		assertEquals(set.getDirection(), DirectionType.AntiSenseDirection);
		assertEquals(set.getNumLists(), 16);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 32);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 1001);

		group = set.getGroup(0);
		assertEquals(group.getDirection(), DirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 16);
		assertEquals(group.getNumCells(), 32);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "Pae_23SrRNA_s_at");

		probe = group.getCell(0);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 124);
		assertEquals(probe.getY(), 95);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("t"), 0);
		probe = group.getCell(1);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 124);
		assertEquals(probe.getY(), 96);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("t"), 0);
		probe = group.getCell(group.getNumCells() - 1);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), group.getNumLists() - 1);
		assertEquals(probe.getX(), 36);
		assertEquals(probe.getY(), 8);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		set = cdf.getProbeSetInformation(header.getNumProbeSets() - 1);
		assertEquals(set.getProbeSetType(), GeneChipProbeSetType.ExpressionProbeSetType);
		assertEquals(set.getDirection(), DirectionType.AntiSenseDirection);
		assertEquals(set.getNumLists(), 20);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 40);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 3101);

		group = set.getGroup(0);
		assertEquals(group.getDirection(), DirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 20);
		assertEquals(group.getNumCells(), 40);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "AFFX_ratb2/X14115_at");

		probe = group.getCell(0);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 20);
		assertEquals(probe.getY(), 113);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		probe = group.getCell(1);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 20);
		assertEquals(probe.getY(), 114);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		probe = group.getCell(group.getNumCells() - 1);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), group.getNumLists() - 1);
		assertEquals(probe.getX(), 119);
		assertEquals(probe.getY(), 58);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		CDFQCProbeSetInformation qcset;
		CDFQCProbeInformation qcprobe;

		qcset = cdf.getQCProbeSetInformation(0);
		assertEquals(qcset.getNumCells(), 300);
		assertEquals(qcset.getQCProbeSetType(), GeneChipQCProbeSetType.GeneExpNegativeQCProbeSetType);

		qcprobe = qcset.getCell(0);
		assertEquals(qcprobe.getX(), 77);
		assertEquals(qcprobe.getY(), 82);
		assertEquals(qcprobe.getProbeLength(), 20);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcprobe = qcset.getCell(1);
		assertEquals(qcprobe.getX(), 77);
		assertEquals(qcprobe.getY(), 83);
		assertEquals(qcprobe.getProbeLength(), 20);
		assertTrue(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcprobe = qcset.getCell(qcset.getNumCells() - 1);
		assertEquals(qcprobe.getX(), 15);
		assertEquals(qcprobe.getY(), 86);
		assertEquals(qcprobe.getProbeLength(), 1);
		assertFalse(qcprobe.isPMProbe());
		assertTrue(qcprobe.isBackground());

		qcset = cdf.getQCProbeSetInformation(header.getNumQCProbeSets() - 1);
		assertEquals(qcset.getNumCells(), 9);
		assertEquals(qcset.getQCProbeSetType(), GeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);

		qcprobe = qcset.getCell(0);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 60);
		assertEquals(qcprobe.getProbeLength(), 25);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcprobe = qcset.getCell(1);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 61);
		assertEquals(qcprobe.getProbeLength(), 25);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcprobe = qcset.getCell(qcset.getNumCells() - 1);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 64);
		assertEquals(qcprobe.getProbeLength(), 25);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcset = cdf.getQCProbeSetInformationByType(GeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);
		assertEquals(qcset.getNumCells(), 9);
		assertEquals(qcset.getQCProbeSetType(), GeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);

		qcprobe = qcset.getCell(0);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 60);
		assertEquals(qcprobe.getProbeLength(), 25);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcprobe = qcset.getCell(1);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 61);
		assertEquals(qcprobe.getProbeLength(), 25);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcprobe = qcset.getCell(qcset.getNumCells() - 1);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 64);
		assertEquals(qcprobe.getProbeLength(), 25);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());
	}

	public void test_ExpressionASCII() {
		CDFFileData cdf = new CDFFileData();
		cdf.setFileName(EXP_ASCII_FILE);
		assertTrue(cdf.read());

		assertEquals(cdf.getChipType(), "Test3-ascii");
		CDFFileHeader header = cdf.getHeader();
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
			assertEquals(cdf.getProbeSetType(i), GeneChipProbeSetType.ExpressionProbeSetType);
		}

		CDFProbeSetInformation set;
		CDFProbeGroupInformation group;
		CDFProbeInformation probe;

		set = cdf.getProbeSetInformation(0);
		assertEquals(set.getProbeSetType(), GeneChipProbeSetType.ExpressionProbeSetType);
		assertEquals(set.getDirection(), DirectionType.AntiSenseDirection);
		assertEquals(set.getNumLists(), 16);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 32);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 1000);

		group = set.getGroup(0);
		assertEquals(group.getDirection(), DirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 16);
		assertEquals(group.getNumCells(), 32);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "Pae_16SrRNA_s_at");

		probe = group.getCell(0);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 111);
		assertEquals(probe.getY(), 79);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		probe = group.getCell(1);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 111);
		assertEquals(probe.getY(), 80);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		probe = group.getCell(group.getNumCells() - 1);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), group.getNumLists() - 1);
		assertEquals(probe.getX(), 93);
		assertEquals(probe.getY(), 82);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("g"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("g"), 0);

		set = cdf.getProbeSetInformation(1);
		assertEquals(set.getProbeSetType(), GeneChipProbeSetType.ExpressionProbeSetType);
		assertEquals(set.getDirection(), DirectionType.AntiSenseDirection);
		assertEquals(set.getNumLists(), 16);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 32);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 1001);

		group = set.getGroup(0);
		assertEquals(group.getDirection(), DirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 16);
		assertEquals(group.getNumCells(), 32);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "Pae_23SrRNA_s_at");

		probe = group.getCell(0);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 124);
		assertEquals(probe.getY(), 95);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("t"), 0);
		probe = group.getCell(1);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 124);
		assertEquals(probe.getY(), 96);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("t"), 0);
		probe = group.getCell(group.getNumCells() - 1);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), group.getNumLists() - 1);
		assertEquals(probe.getX(), 36);
		assertEquals(probe.getY(), 8);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		set = cdf.getProbeSetInformation(header.getNumProbeSets() - 1);
		assertEquals(set.getProbeSetType(), GeneChipProbeSetType.ExpressionProbeSetType);
		assertEquals(set.getDirection(), DirectionType.AntiSenseDirection);
		assertEquals(set.getNumLists(), 20);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 40);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 3101);

		group = set.getGroup(0);
		assertEquals(group.getDirection(), DirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 20);
		assertEquals(group.getNumCells(), 40);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "AFFX_ratb2/X14115_at");

		probe = group.getCell(0);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 20);
		assertEquals(probe.getY(), 113);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		probe = group.getCell(1);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 0);
		assertEquals(probe.getX(), 20);
		assertEquals(probe.getY(), 114);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		probe = group.getCell(group.getNumCells() - 1);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), group.getNumLists() - 1);
		assertEquals(probe.getX(), 119);
		assertEquals(probe.getY(), 58);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		CDFQCProbeSetInformation qcset;
		CDFQCProbeInformation qcprobe;

		qcset = cdf.getQCProbeSetInformation(0);
		assertEquals(qcset.getNumCells(), 300);
		assertEquals(qcset.getQCProbeSetType(), GeneChipQCProbeSetType.GeneExpNegativeQCProbeSetType);

		qcprobe = qcset.getCell(0);
		assertEquals(qcprobe.getX(), 77);
		assertEquals(qcprobe.getY(), 82);
		assertEquals(qcprobe.getProbeLength(), 20);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcprobe = qcset.getCell(1);
		assertEquals(qcprobe.getX(), 77);
		assertEquals(qcprobe.getY(), 83);
		assertEquals(qcprobe.getProbeLength(), 20);

		qcprobe = qcset.getCell(qcset.getNumCells() - 1);
		assertEquals(qcprobe.getX(), 15);
		assertEquals(qcprobe.getY(), 86);
		assertEquals(qcprobe.getProbeLength(), 1);

		qcset = cdf.getQCProbeSetInformation(header.getNumQCProbeSets() - 1);
		assertEquals(qcset.getNumCells(), 9);
		assertEquals(qcset.getQCProbeSetType(), GeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);

		qcprobe = qcset.getCell(0);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 60);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcprobe = qcset.getCell(1);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 61);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcprobe = qcset.getCell(qcset.getNumCells() - 1);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 64);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcset = cdf.getQCProbeSetInformationByType(GeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);
		assertEquals(qcset.getNumCells(), 9);
		assertEquals(qcset.getQCProbeSetType(), GeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);

		qcprobe = qcset.getCell(0);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 60);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcprobe = qcset.getCell(1);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 61);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcprobe = qcset.getCell(qcset.getNumCells() - 1);
		assertEquals(qcprobe.getX(), 62);
		assertEquals(qcprobe.getY(), 64);
		assertEquals(qcprobe.getProbeLength(), 25);
	}

	public void test_GenotypingXDA() {
		CDFFileData cdf = new CDFFileData();
		cdf.setFileName(MAP_XDA_FILE);
		assertTrue(cdf.read());

		assertEquals(cdf.getChipType(), "Mapping10K_Xba131-xda");
		CDFFileHeader header = cdf.getHeader();
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
			assertEquals(cdf.getProbeSetType(i), GeneChipProbeSetType.GenotypingProbeSetType);
		}

		CDFProbeSetInformation set;
		CDFProbeGroupInformation group;
		CDFProbeInformation probe;

		set = cdf.getProbeSetInformation(0);
		assertEquals(set.getProbeSetType(), GeneChipProbeSetType.GenotypingProbeSetType);
		assertEquals(set.getDirection(), DirectionType.SenseDirection);
		assertEquals(set.getNumLists(), 30);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 60);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 41);

		group = set.getGroup(0);
		assertEquals(group.getDirection(), DirectionType.SenseDirection);
		assertEquals(group.getNumLists(), 30);
		assertEquals(group.getNumCells(), 60);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "AFFX-5Q-123");

		probe = group.getCell(0);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 12);
		assertEquals(probe.getX(), 323);
		assertEquals(probe.getY(), 386);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		probe = group.getCell(1);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 12);
		assertEquals(probe.getX(), 323);
		assertEquals(probe.getY(), 387);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		probe = group.getCell(group.getNumCells() - 1);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), 46);
		assertEquals(probe.getX(), 337);
		assertEquals(probe.getY(), 389);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		set = cdf.getProbeSetInformation(header.getNumProbeSets() - 1);
		assertEquals(set.getProbeSetType(), GeneChipProbeSetType.GenotypingProbeSetType);
		assertEquals(set.getDirection(), DirectionType.SenseDirection);
		assertEquals(set.getNumLists(), 20);
		assertEquals(set.getNumGroups(), 4);
		assertEquals(set.getNumCells(), 40);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 12548);

		group = set.getGroup(0);
		assertEquals(group.getDirection(), DirectionType.SenseDirection);
		assertEquals(group.getNumLists(), 5);
		assertEquals(group.getNumCells(), 10);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), 4);
		assertEquals(group.getName(), "SNP_A-1508078A");

		probe = group.getCell(0);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 642);
		assertEquals(probe.getY(), 665);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("g"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("c"), 0);
		probe = group.getCell(1);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 642);
		assertEquals(probe.getY(), 666);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("c"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("c"), 0);
		probe = group.getCell(group.getNumCells() - 1);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), 21);
		assertEquals(probe.getX(), 172);
		assertEquals(probe.getY(), 699);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("t"), 0);

		group = set.getGroup(3);
		assertEquals(group.getDirection(), DirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 5);
		assertEquals(group.getNumCells(), 10);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 15);
		assertEquals(group.getStop(), 19);
		assertEquals(group.getName(), "SNP_A-1508078G");

		probe = group.getCell(0);
		assertEquals(probe.getListIndex(), 15);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 648);
		assertEquals(probe.getY(), 604);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("g"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("g"), 0);
		probe = group.getCell(1);
		assertEquals(probe.getListIndex(), 15);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 648);
		assertEquals(probe.getY(), 603);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("c"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("g"), 0);
		probe = group.getCell(group.getNumCells() - 1);
		assertEquals(probe.getListIndex(), 19);
		assertEquals(probe.getExpos(), 21);
		assertEquals(probe.getX(), 220);
		assertEquals(probe.getY(), 686);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		CDFQCProbeSetInformation qcset;
		CDFQCProbeInformation qcprobe;

		qcset = cdf.getQCProbeSetInformation(0);
		assertEquals(qcset.getNumCells(), 24);
		assertEquals(qcset.getQCProbeSetType(), GeneChipQCProbeSetType.HybPositiveQCProbeSetType);

		qcprobe = qcset.getCell(0);
		assertEquals(qcprobe.getX(), 473);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);
		assertFalse(qcprobe.isPMProbe());
		assertFalse(qcprobe.isBackground());

		qcprobe = qcset.getCell(1);
		assertEquals(qcprobe.getX(), 236);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcprobe = qcset.getCell(qcset.getNumCells() - 1);
		assertEquals(qcprobe.getX(), 711);
		assertEquals(qcprobe.getY(), 472);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcset = cdf.getQCProbeSetInformation(header.getNumQCProbeSets() - 1);
		assertEquals(qcset.getNumCells(), 9);
		assertEquals(qcset.getQCProbeSetType(), GeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);

		qcprobe = qcset.getCell(0);
		assertEquals(qcprobe.getX(), 356);
		assertEquals(qcprobe.getY(), 354);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcprobe = qcset.getCell(1);
		assertEquals(qcprobe.getX(), 356);
		assertEquals(qcprobe.getY(), 355);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcprobe = qcset.getCell(qcset.getNumCells() - 1);
		assertEquals(qcprobe.getX(), 356);
		assertEquals(qcprobe.getY(), 358);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcset = cdf.getQCProbeSetInformationByType(GeneChipQCProbeSetType.HybPositiveQCProbeSetType);
		assertEquals(qcset.getNumCells(), 24);
		assertEquals(qcset.getQCProbeSetType(), GeneChipQCProbeSetType.HybPositiveQCProbeSetType);

		qcprobe = qcset.getCell(0);
		assertEquals(qcprobe.getX(), 473);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcprobe = qcset.getCell(1);
		assertEquals(qcprobe.getX(), 236);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcprobe = qcset.getCell(qcset.getNumCells() - 1);
		assertEquals(qcprobe.getX(), 711);
		assertEquals(qcprobe.getY(), 472);
		assertEquals(qcprobe.getProbeLength(), 16);

	}

	public void test_GenotypingASCII() {
		CDFFileData cdf = new CDFFileData();
		cdf.setFileName(MAP_ASCII_FILE);
		assertTrue(cdf.read());

		assertEquals(cdf.getChipType(), "Mapping10K_Xba131-ascii");
		CDFFileHeader header = cdf.getHeader();
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
			assertEquals(cdf.getProbeSetType(i), GeneChipProbeSetType.GenotypingProbeSetType);
		}

		CDFProbeSetInformation set;
		CDFProbeGroupInformation group;
		CDFProbeInformation probe;

		set = cdf.getProbeSetInformation(0);
		assertEquals(set.getProbeSetType(), GeneChipProbeSetType.GenotypingProbeSetType);
		assertEquals(set.getDirection(), DirectionType.SenseDirection);
		assertEquals(set.getNumLists(), 30);
		assertEquals(set.getNumGroups(), 1);
		assertEquals(set.getNumCells(), 60);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 41);

		group = set.getGroup(0);
		assertEquals(group.getDirection(), DirectionType.SenseDirection);
		assertEquals(group.getNumLists(), 30);
		assertEquals(group.getNumCells(), 60);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), group.getNumLists() - 1);
		assertEquals(group.getName(), "AFFX-5Q-123");

		probe = group.getCell(0);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 12);
		assertEquals(probe.getX(), 323);
		assertEquals(probe.getY(), 386);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("t"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		probe = group.getCell(1);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 12);
		assertEquals(probe.getX(), 323);
		assertEquals(probe.getY(), 387);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);
		probe = group.getCell(group.getNumCells() - 1);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), 46);
		assertEquals(probe.getX(), 337);
		assertEquals(probe.getY(), 389);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		set = cdf.getProbeSetInformation(header.getNumProbeSets() - 1);
		assertEquals(set.getProbeSetType(), GeneChipProbeSetType.GenotypingProbeSetType);
		assertEquals(set.getDirection(), DirectionType.SenseDirection);
		assertEquals(set.getNumLists(), 20);
		assertEquals(set.getNumGroups(), 4);
		assertEquals(set.getNumCells(), 40);
		assertEquals(set.getNumCellsPerList(), 2);
		assertEquals(set.getProbeSetNumber(), 12548);

		group = set.getGroup(0);
		assertEquals(group.getDirection(), DirectionType.SenseDirection);
		assertEquals(group.getNumLists(), 5);
		assertEquals(group.getNumCells(), 10);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 0);
		assertEquals(group.getStop(), 4);
		assertEquals(group.getName(), "SNP_A-1508078A");

		probe = group.getCell(0);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 642);
		assertEquals(probe.getY(), 665);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("g"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("c"), 0);
		probe = group.getCell(1);
		assertEquals(probe.getListIndex(), 0);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 642);
		assertEquals(probe.getY(), 666);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("c"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("c"), 0);
		probe = group.getCell(group.getNumCells() - 1);
		assertEquals(probe.getListIndex(), group.getNumLists() - 1);
		assertEquals(probe.getExpos(), 21);
		assertEquals(probe.getX(), 172);
		assertEquals(probe.getY(), 699);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("t"), 0);

		group = set.getGroup(3);
		assertEquals(group.getDirection(), DirectionType.AntiSenseDirection);
		assertEquals(group.getNumLists(), 5);
		assertEquals(group.getNumCells(), 10);
		assertEquals(group.getNumCellsPerList(), 2);
		assertEquals(group.getStart(), 15);
		assertEquals(group.getStop(), 19);
		assertEquals(group.getName(), "SNP_A-1508078G");

		probe = group.getCell(0);
		assertEquals(probe.getListIndex(), 15);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 648);
		assertEquals(probe.getY(), 604);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("g"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("g"), 0);
		probe = group.getCell(1);
		assertEquals(probe.getListIndex(), 15);
		assertEquals(probe.getExpos(), 15);
		assertEquals(probe.getX(), 648);
		assertEquals(probe.getY(), 603);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("c"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("g"), 0);
		probe = group.getCell(group.getNumCells() - 1);
		assertEquals(probe.getListIndex(), 19);
		assertEquals(probe.getExpos(), 21);
		assertEquals(probe.getX(), 220);
		assertEquals(probe.getY(), 686);
		assertEquals(String.valueOf(probe.getPBase()).compareToIgnoreCase("a"), 0);
		assertEquals(String.valueOf(probe.getTBase()).compareToIgnoreCase("a"), 0);

		CDFQCProbeSetInformation qcset;
		CDFQCProbeInformation qcprobe;

		qcset = cdf.getQCProbeSetInformation(0);
		assertEquals(qcset.getNumCells(), 24);
		assertEquals(qcset.getQCProbeSetType(), GeneChipQCProbeSetType.HybPositiveQCProbeSetType);

		qcprobe = qcset.getCell(0);
		assertEquals(qcprobe.getX(), 473);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcprobe = qcset.getCell(1);
		assertEquals(qcprobe.getX(), 236);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcprobe = qcset.getCell(qcset.getNumCells() - 1);
		assertEquals(qcprobe.getX(), 711);
		assertEquals(qcprobe.getY(), 472);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcset = cdf.getQCProbeSetInformation(header.getNumQCProbeSets() - 1);
		assertEquals(qcset.getNumCells(), 9);
		assertEquals(qcset.getQCProbeSetType(), GeneChipQCProbeSetType.CentralCrossNegativeQCProbeSetType);

		qcprobe = qcset.getCell(0);
		assertEquals(qcprobe.getX(), 356);
		assertEquals(qcprobe.getY(), 354);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcprobe = qcset.getCell(1);
		assertEquals(qcprobe.getX(), 356);
		assertEquals(qcprobe.getY(), 355);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcprobe = qcset.getCell(qcset.getNumCells() - 1);
		assertEquals(qcprobe.getX(), 356);
		assertEquals(qcprobe.getY(), 358);
		assertEquals(qcprobe.getProbeLength(), 25);

		qcset = cdf.getQCProbeSetInformationByType(GeneChipQCProbeSetType.HybPositiveQCProbeSetType);
		assertEquals(qcset.getNumCells(), 24);
		assertEquals(qcset.getQCProbeSetType(), GeneChipQCProbeSetType.HybPositiveQCProbeSetType);

		qcprobe = qcset.getCell(0);
		assertEquals(qcprobe.getX(), 473);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcprobe = qcset.getCell(1);
		assertEquals(qcprobe.getX(), 236);
		assertEquals(qcprobe.getY(), 2);
		assertEquals(qcprobe.getProbeLength(), 16);

		qcprobe = qcset.getCell(qcset.getNumCells() - 1);
		assertEquals(qcprobe.getX(), 711);
		assertEquals(qcprobe.getY(), 472);
		assertEquals(qcprobe.getProbeLength(), 16);
	}

}
