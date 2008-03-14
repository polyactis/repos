/*
 * CHPFileReaderTest.java
 * JUnit based test
 *
 * Created on October 31, 2005, 9:12 AM
 */

package affymetrix.calvin.parsers;

import java.io.File;
import java.io.IOException;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.data.CHPBackgroundZone;
import affymetrix.calvin.data.CHPData;
import affymetrix.calvin.data.CHPExpressionEntry;
import affymetrix.calvin.data.CHPGenotypeEntry;
import affymetrix.calvin.data.CHPReseqCalls;
import affymetrix.calvin.data.CHPReseqEntry;
import affymetrix.calvin.data.CHPReseqForceCall;
import affymetrix.calvin.data.CHPReseqOrigCall;
import affymetrix.calvin.data.CHPUniversalEntry;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;

/**
 * 
 * @author ljevon
 */
public class CHPFileReaderTest extends TestCase {

	public CHPFileReaderTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(CHPFileReaderTest.class);

		return suite;
	}

	/**
	 * Test of getFilename method, of class affymetrix.calvin.parsers.CHPFileReader.
	 */
	public void testFilename() {
		CHPFileReader reader = new CHPFileReader();
		reader.setFilename("file");
		assertEquals(reader.getFilename(), "file");
	}

	public void testReadCHPExpressionFileTest() throws Exception {
		CHPData data = new CHPData();
		CHPFileReader reader = new CHPFileReader();
		String fn = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_expression_file").getCanonicalPath();
		reader.setFilename(fn);
		try {
			reader.read(data);
		}
		catch (Throwable t) {
			assertTrue(false);
		}
		assertEquals(data.getFilename(), fn);
		assertEquals(data.getEntryCount(), 10);
		assertEquals(data.getBackgroundZoneCnt(), 2);

		float eps = 0.000001f;
		CHPExpressionEntry entry1 = new CHPExpressionEntry();
		data.getEntry(0, entry1);
		assertEquals(entry1.getProbeSetName(), "probe set 1");
		assertEquals(entry1.getDetection().toShort(), 10);
		assertEquals(entry1.getDetectionPValue(), 11.0f, eps);
		float signal = entry1.getSignal();
		assertEquals(signal, 17.8f, eps);
		assertEquals(entry1.getNumPairs().toInt(), 6);
		assertEquals(entry1.getNumPairsUsed().toInt(), 5);
		assertEquals(entry1.getChange().toShort(), 2);
		float chgPVal = entry1.getChangePValue();
		assertEquals(chgPVal, 56.9f, eps);
		assertEquals(entry1.getSigLogRatio(), 45.89f, eps);
		assertEquals(entry1.getSigLogRatioLo(), 42.0f, eps);
		assertEquals(entry1.getSigLogRatioHi(), 47.0f, eps);
		assertEquals(entry1.getCommonPairs().toInt(), 2);

		CHPExpressionEntry entry2 = new CHPExpressionEntry();
		data.getEntry(1, entry2);
		assertEquals(entry2.getProbeSetName(), "probe set 2");
		assertEquals(entry2.getDetection().toShort(), 10);
		assertEquals(entry2.getDetectionPValue(), 1.0f, eps);
		assertEquals(entry2.getSignal(), 7.8f, eps);
		assertEquals(entry2.getNumPairs().toInt(), 6);
		assertEquals(entry2.getNumPairsUsed().toInt(), 5);
		assertEquals(entry2.getChange().toShort(), 2);
		assertEquals(entry2.getChangePValue(), 5.9f, eps);
		assertEquals(entry2.getSigLogRatio(), 4.89f, eps);
		assertEquals(entry2.getSigLogRatioLo(), 2.0f, eps);
		assertEquals(entry2.getSigLogRatioHi(), 7.0f, eps);
		assertEquals(entry2.getCommonPairs().toInt(), 2);

		CHPBackgroundZone zone1 = new CHPBackgroundZone();
		data.getBackgroundZone(0, zone1);
		assertEquals(zone1.getCenterX(), 11.0f, eps);
		assertEquals(zone1.getCenterY(), 17.8f, eps);
		assertEquals(zone1.getBackground(), 56.9f, eps);
		assertEquals(zone1.getSmoothFactor(), 45.89f, eps);

		CHPBackgroundZone zone2 = new CHPBackgroundZone();
		data.getBackgroundZone(1, zone2);
		assertEquals(zone2.getCenterX(), 1.0f, eps);
		assertEquals(zone2.getCenterY(), 7.8f, eps);
		assertEquals(zone2.getBackground(), 6.9f, eps);
		assertEquals(zone2.getSmoothFactor(), 5.89f, eps);
	}

	public void testReadCHPGenotypingFileTest() throws IOException, UnsignedOutOfLimitsException {
		CHPData data = new CHPData();
		CHPFileReader reader = new CHPFileReader();
		String fn = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_genotype_file").getCanonicalPath();
		reader.setFilename(fn);
		try {
			reader.read(data);
		}
		catch (Throwable t) {
			assertTrue(false);
		}

		assertEquals(data.getFilename(), fn);
		assertEquals(data.getEntryCount(), 2);
		assertEquals(data.getBackgroundZoneCnt(), 2);

		float eps = 0.00001f;
		CHPGenotypeEntry entry1 = new CHPGenotypeEntry();
		data.getEntry(0, entry1);
		assertEquals(entry1.getProbeSetName(), "probe set 1");
		assertEquals(entry1.getCall(), 1);
		assertEquals(entry1.getConfidence(), 11.0f, eps);
		float ras1 = entry1.getRAS1();
		assertEquals(ras1, 17.8f, eps);
		assertEquals(entry1.getRAS2(), 6.0f, eps);
		assertEquals(entry1.getAACall(), 5.0f, eps);
		assertEquals(entry1.getABCall(), 2.0f, eps);
		assertEquals(entry1.getBBCall(), 56.9f, eps);
		assertEquals(entry1.getNoCall(), 45.89f, eps);

		CHPGenotypeEntry entry2 = new CHPGenotypeEntry();
		data.getEntry(1, entry2);
		assertEquals(entry2.getProbeSetName(), "probe set 2");
		assertEquals(entry2.getCall(), 2);
		assertEquals(entry2.getConfidence(), 1.0f, eps);
		assertEquals(entry2.getRAS1(), 7.8f, eps);
		assertEquals(entry2.getRAS2(), 6.0f, eps);
		assertEquals(entry2.getAACall(), 5.0f, eps);
		assertEquals(entry2.getABCall(), 2.0f, eps);
		assertEquals(entry2.getBBCall(), 5.9f, eps);
		assertEquals(entry2.getNoCall(), 4.8f, eps);

		CHPBackgroundZone zone1 = new CHPBackgroundZone();
		data.getBackgroundZone(0, zone1);
		assertEquals(zone1.getCenterX(), 11.0f, eps);
		assertEquals(zone1.getCenterY(), 17.8f, eps);
		assertEquals(zone1.getBackground(), 56.9f, eps);
		assertEquals(zone1.getSmoothFactor(), 45.89f, eps);

		CHPBackgroundZone zone2 = new CHPBackgroundZone();
		data.getBackgroundZone(1, zone2);
		assertEquals(zone2.getCenterX(), 1.0f, eps);
		assertEquals(zone2.getCenterY(), 7.8f, eps);
		assertEquals(zone2.getBackground(), 6.9f, eps);
		assertEquals(zone2.getSmoothFactor(), 5.89f, eps);
	}

	public void testReadCHPUniversalFileTest() throws IOException, UnsignedOutOfLimitsException {
		CHPData data = new CHPData();
		CHPFileReader reader = new CHPFileReader();
		String fn = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_universal_file").getCanonicalPath();
		reader.setFilename(fn);
		try {
			reader.read(data);
		}
		catch (Throwable t) {
			assertTrue(false);
		}

		assertEquals(data.getFilename(), fn);
		assertEquals(data.getEntryCount(), 2);
		assertEquals(data.getBackgroundZoneCnt(), 2);

		float eps = 0.00001f;
		CHPUniversalEntry entry1 = new CHPUniversalEntry();
		data.getEntry(0, entry1);
		assertEquals(entry1.getBackground(), 11.0f, eps);

		CHPUniversalEntry entry2 = new CHPUniversalEntry();
		data.getEntry(1, entry2);
		assertEquals(entry2.getBackground(), 1.0f, eps);

		CHPBackgroundZone zone1 = new CHPBackgroundZone();
		data.getBackgroundZone(0, zone1);
		assertEquals(zone1.getCenterX(), 11.0f, eps);
		assertEquals(zone1.getCenterY(), 17.8f, eps);
		assertEquals(zone1.getBackground(), 56.9f, eps);
		assertEquals(zone1.getSmoothFactor(), 45.89f, eps);

		CHPBackgroundZone zone2 = new CHPBackgroundZone();
		data.getBackgroundZone(1, zone2);
		assertEquals(zone2.getCenterX(), 1.0f, eps);
		assertEquals(zone2.getCenterY(), 7.8f, eps);
		assertEquals(zone2.getBackground(), 6.9f, eps);
		assertEquals(zone2.getSmoothFactor(), 5.89f, eps);
	}

	public void testReadCHPReseqFileTest() throws IOException, UnsignedOutOfLimitsException {
		CHPData data = new CHPData();
		CHPFileReader reader = new CHPFileReader();
		String fn = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_reseq_file").getCanonicalPath();
		reader.setFilename(fn);
		try {
			reader.read(data);
		}
		catch (Throwable t) {
			assertTrue(false);
		}

		assertEquals(data.getFilename(), fn);

		assertEquals(data.getEntryCount(), 5);
		assertEquals(data.getBackgroundZoneCnt(), 1);
		assertEquals(data.getForceCnt(), 2);
		assertEquals(data.getOrigCnt(), 3);

		float eps = 0.000001f;
		CHPReseqEntry e = new CHPReseqEntry();
		data.getEntry(0, e);
		assertEquals(e.getCall(), (byte)'a');
		assertEquals(e.getScore(), 1.0f, eps);

		data.getEntry(1, e);
		assertEquals(e.getCall(), (byte)'c');
		assertEquals(e.getScore(), 2.0f, eps);

		data.getEntry(2, e);
		assertEquals(e.getCall(), (byte)'g');
		assertEquals(e.getScore(), 3.0f, eps);

		data.getEntry(3, e);
		assertEquals(e.getCall(), (byte)'t');
		assertEquals(e.getScore(), 4.0f, eps);

		data.getEntry(4, e);
		assertEquals(e.getCall(), (byte)'n');
		assertEquals(e.getScore(), 5.0f, eps);

		CHPReseqForceCall force = new CHPReseqForceCall();
		data.getForceCall(0, force);
		assertEquals(force.getPosition(), 1);
		assertEquals(force.getCall(), (byte)'a');
		assertEquals(force.getReason(), CHPReseqCalls.CC_SATURATION_LEVEL_FORCE_CALL);

		data.getForceCall(1, force);
		assertEquals(force.getPosition(), 2);
		assertEquals(force.getCall(), (byte)'c');
		assertEquals(force.getReason(), CHPReseqCalls.CC_WEAK_SIGNAL_THR_FORCE_CALL);

		CHPReseqOrigCall orig = new CHPReseqOrigCall();
		data.getOrigCall(0, orig);
		assertEquals(orig.getPosition(), 3);
		assertEquals(orig.getCall(), (byte)'t');

		data.getOrigCall(1, orig);
		assertEquals(orig.getPosition(), 4);
		assertEquals(orig.getCall(), (byte)'a');

		data.getOrigCall(2, orig);
		assertEquals(orig.getPosition(), 5);
		assertEquals(orig.getCall(), (byte)'g');

	}

}