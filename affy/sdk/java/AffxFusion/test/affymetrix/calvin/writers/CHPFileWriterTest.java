package affymetrix.calvin.writers;

import junit.framework.TestCase;
import affymetrix.calvin.data.CHPBackgroundZone;
import affymetrix.calvin.data.CHPData;
import affymetrix.calvin.data.CHPExpressionEntry;
import affymetrix.calvin.data.CHPGenotypeEntry;
import affymetrix.calvin.data.CHPReseqCalls;
import affymetrix.calvin.data.CHPReseqEntry;
import affymetrix.calvin.data.CHPReseqForceCall;
import affymetrix.calvin.data.CHPReseqOrigCall;
import affymetrix.calvin.data.CHPUniversalEntry;
import affymetrix.calvin.parsers.CHPFileReader;
import affymetrix.portability.UByte;
import affymetrix.portability.UShort;

public class CHPFileWriterTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testWriteExpressionEntry() throws Exception {
		String filename = "CHP_expression_file";
		CHPData data = new CHPData(filename, CHPData.CHP_EXPRESSION_ASSAY_TYPE);
		data.setEntryCount(2, 16, true);
		data.setBackgroundZoneCnt(2);
		CHPFileWriter writer = new CHPFileWriter(data);

		CHPExpressionEntry e1 = new CHPExpressionEntry("probe set 1", new UByte((byte)10), 11.0f, 17.8f, new UShort(6),
				new UShort(5), true, new UByte((byte)2), 56.9f, 45.89f, 42.0f, 47.0f, new UShort(2));
		CHPExpressionEntry e2 = new CHPExpressionEntry("probe set 2", new UByte((byte)10), 1.0f, 7.8f, new UShort(6),
				new UShort(5), true, new UByte((byte)2), 5.9f, 4.89f, 2.0f, 7.0f, new UShort(2));
		writer.seekToDataSet();
		writer.writeExpressionEntry(e1);
		writer.writeExpressionEntry(e2);

		CHPBackgroundZone bgz1 = new CHPBackgroundZone(11.0f, 17.8f, 56.9f, 45.89f);
		CHPBackgroundZone bgz2 = new CHPBackgroundZone(1.0f, 7.8f, 6.9f, 5.89f);
		writer.seekToBgSet();
		writer.writeBackgroundZone(bgz1);
		writer.writeBackgroundZone(bgz2);
		writer.close();
		assertTrue(true);

		CHPFileReader reader = new CHPFileReader();
		reader.setFilename(filename);
		CHPData data2 = new CHPData();
		reader.read(data2);

		CHPExpressionEntry e = new CHPExpressionEntry();
		data2.getEntry(0, e);
		assertEquals(e.getProbeSetName(), e1.getProbeSetName());
		assertEquals(e.getSignal(), e1.getSignal(), 0.000001f);
		assertEquals(e.getDetectionPValue(), e1.getDetectionPValue(), 0.000001f);
		assertEquals(e.getChangePValue(), e1.getChangePValue(), 0.000001f);
		assertEquals(e.getSigLogRatio(), e1.getSigLogRatio(), 0.000001f);
		assertEquals(e.getSigLogRatioLo(), e1.getSigLogRatioLo(), 0.000001f);
		assertEquals(e.getSigLogRatioHi(), e1.getSigLogRatioHi(), 0.000001f);
		assertEquals(e.getDetection().toShort(), e1.getDetection().toShort());
		assertEquals(e.getNumPairs().toInt(), e1.getNumPairs().toInt());
		assertEquals(e.getNumPairsUsed().toInt(), e1.getNumPairsUsed().toInt());
		assertEquals(e.getHasComparisonData(), e1.getHasComparisonData());
		assertEquals(e.getChange().toShort(), e1.getChange().toShort());
		assertEquals(e.getCommonPairs().toInt(), e1.getCommonPairs().toInt());

		data2.getEntry(1, e);
		assertEquals(e.getProbeSetName(), e2.getProbeSetName());
		assertEquals(e.getSignal(), e2.getSignal(), 0.000001f);
		assertEquals(e.getDetectionPValue(), e2.getDetectionPValue(), 0.000001f);
		assertEquals(e.getChangePValue(), e2.getChangePValue(), 0.000001f);
		assertEquals(e.getSigLogRatio(), e2.getSigLogRatio(), 0.000001f);
		assertEquals(e.getSigLogRatioLo(), e2.getSigLogRatioLo(), 0.000001f);
		assertEquals(e.getSigLogRatioHi(), e2.getSigLogRatioHi(), 0.000001f);
		assertEquals(e.getDetection().toShort(), e2.getDetection().toShort());
		assertEquals(e.getNumPairs().toInt(), e2.getNumPairs().toInt());
		assertEquals(e.getNumPairsUsed().toInt(), e2.getNumPairsUsed().toInt());
		assertEquals(e.getHasComparisonData(), e2.getHasComparisonData());
		assertEquals(e.getChange().toShort(), e2.getChange().toShort());
		assertEquals(e.getCommonPairs().toInt(), e2.getCommonPairs().toInt());

		CHPBackgroundZone bgz = new CHPBackgroundZone();
		data2.getBackgroundZone(0, bgz);
		assertEquals(bgz.getBackground(), bgz1.getBackground(), 0.000001f);
		assertEquals(bgz.getCenterX(), bgz1.getCenterX(), 0.000001f);
		assertEquals(bgz.getCenterY(), bgz1.getCenterY(), 0.000001f);
		assertEquals(bgz.getSmoothFactor(), bgz1.getSmoothFactor(), 0.000001f);

		data2.getBackgroundZone(1, bgz);
		assertEquals(bgz.getBackground(), bgz2.getBackground(), 0.000001f);
		assertEquals(bgz.getCenterX(), bgz2.getCenterX(), 0.000001f);
		assertEquals(bgz.getCenterY(), bgz2.getCenterY(), 0.000001f);
		assertEquals(bgz.getSmoothFactor(), bgz2.getSmoothFactor(), 0.000001f);
	}

	public void testWriteGenotypeEntry() throws Exception {
		CHPData data = new CHPData("CHP_genotype_file", CHPData.CHP_GENOTYPING_ASSAY_TYPE);
		data.setEntryCount(2, 16, false);
		data.setBackgroundZoneCnt(2);
		CHPFileWriter writer = new CHPFileWriter(data);

		CHPGenotypeEntry e1 = new CHPGenotypeEntry("probe set 1", (byte)1, 11.0f, 17.8f, 6.0f, 5.0f, 2.0f, 56.9f, 45.89f);
		CHPGenotypeEntry e2 = new CHPGenotypeEntry("probe set 2", (byte)2, 1.0f, 7.8f, 6.0f, 5.0f, 2.0f, 5.9f, 4.8f);
		writer.seekToDataSet();
		writer.writeGenotypeEntry(e1);
		writer.writeGenotypeEntry(e2);

		CHPBackgroundZone bgz1 = new CHPBackgroundZone(11.0f, 17.8f, 56.9f, 45.89f);
		CHPBackgroundZone bgz2 = new CHPBackgroundZone(1.0f, 7.8f, 6.9f, 5.89f);
		writer.seekToBgSet();
		writer.writeBackgroundZone(bgz1);
		writer.writeBackgroundZone(bgz2);
		writer.close();
		assertTrue(true);

		CHPFileReader reader = new CHPFileReader();
		reader.setFilename("CHP_genotype_file");
		CHPData data2 = new CHPData();
		reader.read(data2);

		CHPGenotypeEntry e = new CHPGenotypeEntry();
		data2.getEntry(0, e);
		assertEquals(e.getProbeSetName(), e1.getProbeSetName());
		assertTrue(e.getCall() == e1.getCall());
		assertEquals(e.getConfidence(), e1.getConfidence(), 0.00001f);
		assertEquals(e.getRAS1(), e1.getRAS1(), 0.00001f);
		assertEquals(e.getRAS2(), e1.getRAS2(), 0.00001f);
		assertEquals(e.getAACall(), e1.getAACall(), 0.00001f);
		assertEquals(e.getABCall(), e1.getABCall(), 0.00001f);
		assertEquals(e.getBBCall(), e1.getBBCall(), 0.00001f);
		assertEquals(e.getNoCall(), e1.getNoCall(), 0.00001f);

		data2.getEntry(1, e);
		assertEquals(e.getProbeSetName(), e2.getProbeSetName());
		assertTrue(e.getCall() == e2.getCall());
		assertEquals(e.getConfidence(), e2.getConfidence(), 0.00001f);
		assertEquals(e.getRAS1(), e2.getRAS1(), 0.00001f);
		assertEquals(e.getRAS2(), e2.getRAS2(), 0.00001f);
		assertEquals(e.getAACall(), e2.getAACall(), 0.00001f);
		assertEquals(e.getABCall(), e2.getABCall(), 0.00001f);
		assertEquals(e.getBBCall(), e2.getBBCall(), 0.00001f);
		assertEquals(e.getNoCall(), e2.getNoCall(), 0.00001f);

		CHPBackgroundZone bgz = new CHPBackgroundZone();
		data2.getBackgroundZone(0, bgz);
		assertEquals(bgz.getBackground(), bgz1.getBackground(), 0.000001f);
		assertEquals(bgz.getCenterX(), bgz1.getCenterX(), 0.000001f);
		assertEquals(bgz.getCenterY(), bgz1.getCenterY(), 0.000001f);
		assertEquals(bgz.getSmoothFactor(), bgz1.getSmoothFactor(), 0.000001f);

		data2.getBackgroundZone(1, bgz);
		assertEquals(bgz.getBackground(), bgz2.getBackground(), 0.000001f);
		assertEquals(bgz.getCenterX(), bgz2.getCenterX(), 0.000001f);
		assertEquals(bgz.getCenterY(), bgz2.getCenterY(), 0.000001f);
		assertEquals(bgz.getSmoothFactor(), bgz2.getSmoothFactor(), 0.000001f);

	}

	public void testWriteUniversalEntry() throws Exception {
		CHPData data = new CHPData("CHP_universal_file", CHPData.CHP_UNIVERSAL_ASSAY_TYPE);
		data.setEntryCount(2, 8, false);
		data.setBackgroundZoneCnt(2);
		CHPFileWriter writer = new CHPFileWriter(data);

		CHPUniversalEntry e1 = new CHPUniversalEntry(11.0f);
		CHPUniversalEntry e2 = new CHPUniversalEntry(1.0f);
		writer.seekToDataSet();
		writer.writeUniversalEntry(e1);
		writer.writeUniversalEntry(e2);

		CHPBackgroundZone bgz1 = new CHPBackgroundZone(11.0f, 17.8f, 56.9f, 45.89f);
		CHPBackgroundZone bgz2 = new CHPBackgroundZone(1.0f, 7.8f, 6.9f, 5.89f);
		writer.seekToBgSet();
		writer.writeBackgroundZone(bgz1);
		writer.writeBackgroundZone(bgz2);
		writer.close();
		assertTrue(true);

		CHPFileReader reader = new CHPFileReader();
		reader.setFilename("CHP_universal_file");
		CHPData data2 = new CHPData();
		reader.read(data2);

		CHPUniversalEntry e = new CHPUniversalEntry();

		data2.getEntry(0, e);
		assertEquals(e.getBackground(), e1.getBackground(), 0.000001f);
		data2.getEntry(1, e);
		assertEquals(e.getBackground(), e2.getBackground(), 0.000001f);

		CHPBackgroundZone bgz = new CHPBackgroundZone();
		data2.getBackgroundZone(0, bgz);
		assertEquals(bgz.getBackground(), bgz1.getBackground(), 0.000001f);
		assertEquals(bgz.getCenterX(), bgz1.getCenterX(), 0.000001f);
		assertEquals(bgz.getCenterY(), bgz1.getCenterY(), 0.000001f);
		assertEquals(bgz.getSmoothFactor(), bgz1.getSmoothFactor(), 0.000001f);

		data2.getBackgroundZone(1, bgz);
		assertEquals(bgz.getBackground(), bgz2.getBackground(), 0.000001f);
		assertEquals(bgz.getCenterX(), bgz2.getCenterX(), 0.000001f);
		assertEquals(bgz.getCenterY(), bgz2.getCenterY(), 0.000001f);
		assertEquals(bgz.getSmoothFactor(), bgz2.getSmoothFactor(), 0.000001f);
	}

	public void testWriteReseqEntry() throws Exception {
		CHPData data = new CHPData("CHP_reseq_file", CHPData.CHP_RESEQUENCING_ASSAY_TYPE);
		data.setEntryCount(5, 10, false);
		data.setBackgroundZoneCnt(1);
		data.setForceCnt(2);
		data.setOrigCnt(3);
		CHPFileWriter writer = new CHPFileWriter(data);
		writer.seekToDataSet();

		CHPReseqEntry e = new CHPReseqEntry();
		e.setCall((byte)'a');
		e.setScore(1.0f);
		writer.writeReseqEntry(e);
		e.setCall((byte)'c');
		e.setScore(2.0f);
		writer.writeReseqEntry(e);
		e.setCall((byte)'g');
		e.setScore(3.0f);
		writer.writeReseqEntry(e);
		e.setCall((byte)'t');
		e.setScore(4.0f);
		writer.writeReseqEntry(e);
		e.setCall((byte)'n');
		e.setScore(5.0f);
		writer.writeReseqEntry(e);

		writer.seekToForceSet();
		CHPReseqForceCall force = new CHPReseqForceCall();
		force.setPosition(1);
		force.setCall((byte)'a');
		force.setReason(CHPReseqCalls.CC_SATURATION_LEVEL_FORCE_CALL);
		writer.writeForceCall(force);
		force.setPosition(2);
		force.setCall((byte)'c');
		force.setReason(CHPReseqCalls.CC_WEAK_SIGNAL_THR_FORCE_CALL);
		writer.writeForceCall(force);

		writer.seekToOrigCallSet();
		CHPReseqOrigCall orig = new CHPReseqOrigCall();
		orig.setPosition(3);
		orig.setCall((byte)'t');
		writer.writeOrigCall(orig);
		orig.setPosition(4);
		orig.setCall((byte)'a');
		writer.writeOrigCall(orig);
		orig.setPosition(5);
		orig.setCall((byte)'g');
		writer.writeOrigCall(orig);
		writer.close();
		assertTrue(true);

		CHPFileReader reader = new CHPFileReader();
		reader.setFilename("CHP_reseq_file");
		CHPData data2 = new CHPData();
		reader.read(data2);

		data2.getEntry(0, e);
		assertTrue(e.getCall() == 'a');
		assertEquals(e.getScore(), 1.0f, 0.0001f);

		data2.getEntry(1, e);
		assertTrue(e.getCall() == 'c');
		assertEquals(e.getScore(), 2.0f, 0.0001f);

		data2.getEntry(2, e);
		assertTrue(e.getCall() == 'g');
		assertEquals(e.getScore(), 3.0f, 0.0001f);

		data2.getEntry(3, e);
		assertTrue(e.getCall() == 't');
		assertEquals(e.getScore(), 4.0f, 0.0001f);

		data2.getEntry(4, e);
		assertTrue(e.getCall() == 'n');
		assertEquals(e.getScore(), 5.0f, 0.0001f);

		data2.getForceCall(0, force);
		assertTrue(force.getCall() == 'a');
		assertTrue(force.getReason() == CHPReseqCalls.CC_SATURATION_LEVEL_FORCE_CALL);
		assertTrue(force.getPosition() == 1);

		data2.getForceCall(1, force);
		assertTrue(force.getCall() == 'c');
		assertTrue(force.getReason() == CHPReseqCalls.CC_WEAK_SIGNAL_THR_FORCE_CALL);
		assertTrue(force.getPosition() == 2);

		data2.getOrigCall(0, orig);
		assertTrue(orig.getCall() == 't');
		assertTrue(orig.getPosition() == 3);

		data2.getOrigCall(1, orig);
		assertTrue(orig.getCall() == 'a');
		assertTrue(orig.getPosition() == 4);

		data2.getOrigCall(2, orig);
		assertTrue(orig.getCall() == 'g');
		assertTrue(orig.getPosition() == 5);
	}
}
