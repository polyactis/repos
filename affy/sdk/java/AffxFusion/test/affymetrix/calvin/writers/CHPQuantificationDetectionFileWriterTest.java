package affymetrix.calvin.writers;

import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;
import affymetrix.calvin.data.CHPQuantificationDetectionData;
import affymetrix.calvin.data.ProbeSetQuantificationDetectionData;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.parsers.CHPQuantificationDetectionFileReader;

public class CHPQuantificationDetectionFileWriterTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testWrite() throws Exception {
		List<ParameterNameValue> params = new ArrayList<ParameterNameValue>();

		CHPQuantificationDetectionData data = new CHPQuantificationDetectionData("CHP_quantification_detection_file");

		data.setAlgName("sig");
		data.setAlgVersion("1.0");
		data.setArrayType("test3");
		data.setEntryCount(2, 10);

		ParameterNameValue param1 = new ParameterNameValue();
		param1.setName("an1");
		param1.setValueText("av1");
		params.add(param1);
		data.addAlgParams(params);

		params.clear();
		ParameterNameValue param2 = new ParameterNameValue();
		param2.setName("sn1");
		param2.setValueText("sv1");
		params.add(param2);
		data.addSummaryParams(params);

		CHPQuantificationDetectionFileWriter writer = new CHPQuantificationDetectionFileWriter(data);
		ProbeSetQuantificationDetectionData e = new ProbeSetQuantificationDetectionData();
		writer.seekToDataSet();
		e.setName("abc");
		e.setQuantification(10f);
		e.setPValue(.1f);
		writer.writeEntry(e);
		e.setName("xyz");
		e.setQuantification(20f);
		e.setPValue(0.2f);
		writer.writeEntry(e);
		writer.close();
		assertTrue(true);

		CHPQuantificationDetectionData data2 = new CHPQuantificationDetectionData();
		CHPQuantificationDetectionFileReader reader = new CHPQuantificationDetectionFileReader();
		reader.setFilename("CHP_quantification_detection_file");
		reader.read(data2);

		assertEquals(data2.getAlgName(), "sig");
		assertEquals(data2.getAlgVersion(), "1.0");
		assertEquals(data2.getArrayType(), "test3");
		assertEquals(data2.getEntryCount(), 2);

		ParameterNameValue param = null;
		param = data2.getAlgParams().iterator().next();
		assertEquals(param.getName(), "an1");
		assertEquals(param.getValueText(), "av1");

		param = data2.getSummaryParams().iterator().next();
		assertEquals(param.getName(), "sn1");
		assertEquals(param.getValueText(), "sv1");

		e = data2.getQuantificationDetectionEntry(0);
		assertEquals(e.getQuantification(), 10.0f, 0.0001f);
		assertEquals(e.getPValue(), 0.10f, 0.0001f);
		assertEquals(e.getName(), "abc");
		e = data2.getQuantificationDetectionEntry(1);
		assertEquals(e.getQuantification(), 20.0f, 0.0001f);
		assertEquals(e.getPValue(), 0.20f, 0.0001f);
		assertEquals(e.getName(), "xyz");

	}

	public void testWriteId() throws Exception {
		List<ParameterNameValue> params = new ArrayList<ParameterNameValue>();

		CHPQuantificationDetectionData data = new CHPQuantificationDetectionData("CHP_quantification_detection_file_id");

		data.setAlgName("sig");
		data.setAlgVersion("1.0");
		data.setArrayType("test3");
		data.setEntryCount(2);

		ParameterNameValue param1 = new ParameterNameValue();
		param1.setName("an1");
		param1.setValueText("av1");
		params.add(param1);
		data.addAlgParams(params);

		params.clear();
		ParameterNameValue param2 = new ParameterNameValue();
		param2.setName("sn1");
		param2.setValueText("sv1");
		params.add(param2);
		data.addSummaryParams(params);

		CHPQuantificationDetectionFileWriter writer = new CHPQuantificationDetectionFileWriter(data);
		ProbeSetQuantificationDetectionData e = new ProbeSetQuantificationDetectionData();
		writer.seekToDataSet();
		e.setId(10);
		e.setQuantification(10f);
		e.setPValue(.1f);
		writer.writeEntry(e);
		e.setId(20);
		e.setQuantification(20f);
		e.setPValue(.2f);
		writer.writeEntry(e);
		writer.close();
		assertTrue(true);

		CHPQuantificationDetectionData data2 = new CHPQuantificationDetectionData();
		CHPQuantificationDetectionFileReader reader = new CHPQuantificationDetectionFileReader();
		reader.setFilename("CHP_quantification_detection_file_id");
		reader.read(data2);

		assertEquals(data2.getAlgName(), "sig");
		assertEquals(data2.getAlgVersion(), "1.0");
		assertEquals(data2.getArrayType(), "test3");
		assertEquals(data2.getEntryCount(), 2);

		ParameterNameValue param = null;
		param = data2.getAlgParams().iterator().next();
		assertEquals(param.getName(), "an1");
		assertEquals(param.getValueText(), "av1");

		param = data2.getSummaryParams().iterator().next();
		assertEquals(param.getName(), "sn1");
		assertEquals(param.getValueText(), "sv1");

		e = data2.getQuantificationDetectionEntry(0);
		assertEquals(e.getQuantification(), 10.0f, 0.0001f);
		assertEquals(e.getPValue(), 0.1f, 0.0001f);
		assertEquals(e.getId(), 10);
		assertEquals(e.getName(), "");
		e = data2.getQuantificationDetectionEntry(1);
		assertEquals(e.getQuantification(), 20.0f, 0.0001f);
		assertEquals(e.getPValue(), 0.2f, 0.0001f);
		assertEquals(e.getId(), 20);
		assertEquals(e.getName(), "");

	}
}
