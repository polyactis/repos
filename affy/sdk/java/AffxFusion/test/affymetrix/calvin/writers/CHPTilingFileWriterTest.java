package affymetrix.calvin.writers;

import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;
import affymetrix.calvin.data.CHPTilingData;
import affymetrix.calvin.data.CHPTilingEntry;
import affymetrix.calvin.data.TilingSequenceData;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.parsers.CHPTilingFileReader;

public class CHPTilingFileWriterTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testWrite() throws Exception {
		List<ParameterNameValue> params = new ArrayList<ParameterNameValue>();
		CHPTilingData data = new CHPTilingData("CHP_tiling_file");

		data.setAlgName("tile");
		data.setAlgVersion("1.0");
		data.setNumberSequences(2);
		ParameterNameValue param1 = new ParameterNameValue();
		param1.setName("p1");
		param1.setValueText("v1");
		params.add(param1);
		data.addAlgParams(params);

		TilingSequenceData seq1 = new TilingSequenceData();
		seq1.setName("n1");
		seq1.setGroupName("g1");
		seq1.setVersion("v1");
		ParameterNameValue param2 = new ParameterNameValue();
		param2.setName("seq1_p1");
		param2.setValueText("seq1_v1");
		seq1.addParameter(param2);
		data.addTilingSequenceData(2, seq1);

		TilingSequenceData seq2 = new TilingSequenceData();
		seq2.setName("n2");
		seq2.setGroupName("g2");
		seq2.setVersion("v2");
		ParameterNameValue param3 = new ParameterNameValue();
		param3.setName("seq2_p1");
		param3.setValueText("seq2_v1");
		seq2.addParameter(param3);
		data.addTilingSequenceData(3, seq2);

		CHPTilingFileWriter writer = new CHPTilingFileWriter(data);
		CHPTilingEntry e = new CHPTilingEntry();

		writer.seekToDataSet(0);
		e.setPosition(10);
		e.setValue(10.0f);
		writer.writeTilingEntry(e);
		e.setPosition(20);
		e.setValue(20.0f);
		writer.writeTilingEntry(e);

		writer.seekToDataSet(1);
		e.setPosition(11);
		e.setValue(11.0f);
		writer.writeTilingEntry(e);
		e.setPosition(21);
		e.setValue(21.0f);
		writer.writeTilingEntry(e);
		e.setPosition(31);
		e.setValue(31.0f);
		writer.writeTilingEntry(e);
		writer.close();
		assertTrue(true);

		CHPTilingData data2 = new CHPTilingData();
		CHPTilingFileReader reader = new CHPTilingFileReader();
		reader.setFilename("CHP_tiling_file");
		reader.read(data2);

		assertEquals(data2.getAlgName(), "tile");
		assertEquals(data2.getAlgVersion(), "1.0");
		assertEquals(data2.getNumberSequences(), 2);

		ParameterNameValue param = null;
		param = data2.getAlgParams().iterator().next();
		assertEquals(param.getName(), "p1");
		assertEquals(param.getValueText(), "v1");
		assertEquals(data2.getNumberSequences(), 2);

		data2.openTilingSequenceDataSet(0);
		TilingSequenceData seq = data2.getTilingSequenceData();
		assertEquals(seq.getName(), "n1");
		assertEquals(seq.getGroupName(), "g1");
		assertEquals(seq.getVersion(), "v1");
		param = seq.getParameters().iterator().next();
		assertEquals(param.getName(), "seq1_p1");
		assertEquals(param.getValueText(), "seq1_v1");
		assertEquals(data2.getTilingSequenceEntryCount(0), 2);

		e = data2.getTilingSequenceEntry(0);
		assertEquals(e.getPosition(), 10);
		assertEquals(e.getValue(), 10.0f, 0.00001f);
		e = data2.getTilingSequenceEntry(1);
		assertEquals(e.getPosition(), 20);
		assertEquals(e.getValue(), 20.0f, 0.00001f);

		data2.openTilingSequenceDataSet(1);
		seq = data2.getTilingSequenceData();
		assertEquals(seq.getName(), "n2");
		assertEquals(seq.getGroupName(), "g2");
		assertEquals(seq.getVersion(), "v2");
		param = seq.getParameters().iterator().next();
		assertEquals(param.getName(), "seq2_p1");
		assertEquals(param.getValueText(), "seq2_v1");
		assertEquals(data2.getTilingSequenceEntryCount(1), 3);

		e = data2.getTilingSequenceEntry(0);
		assertEquals(e.getPosition(), 11);
		assertEquals(e.getValue(), 11.0f, 0.00001f);
		e = data2.getTilingSequenceEntry(1);
		assertEquals(e.getPosition(), 21);
		assertEquals(e.getValue(), 21.0f, 0.00001f);
		e = data2.getTilingSequenceEntry(2);
		assertEquals(e.getPosition(), 31);
		assertEquals(e.getValue(), 31.0f, 0.00001f);

	}
}
