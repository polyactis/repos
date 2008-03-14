package affymetrix.calvin.writers;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import junit.framework.TestCase;
import affymetrix.calvin.data.CHPQuantificationData;
import affymetrix.calvin.data.ProbeSetQuantificationData;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.parsers.CHPQuantificationFileReader;

public class CHPQuantificationFileWriterTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testWrite() throws Exception {
		List<ParameterNameValue> params = new ArrayList<ParameterNameValue>();
		CHPQuantificationData data = new CHPQuantificationData("CHP_quantification_file");
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

		CHPQuantificationFileWriter writer = new CHPQuantificationFileWriter(data);
		ProbeSetQuantificationData psQuantData = new ProbeSetQuantificationData();
		writer.seekToDataSet();
		psQuantData.setName("abc");
		psQuantData.setQuantification(10.0f);
		writer.writeEntry(psQuantData);
		psQuantData.setName("xyz");
		psQuantData.setQuantification(20.0f);
		writer.writeEntry(psQuantData);
		writer.close();
		assertTrue(true);

		CHPQuantificationData data2 = new CHPQuantificationData();
		CHPQuantificationFileReader reader = new CHPQuantificationFileReader();
		reader.setFilename("CHP_quantification_file");
		reader.read(data2);

		assertEquals(data2.getAlgName(), "sig");
		assertEquals(data2.getAlgVersion(), "1.0");
		assertEquals(data2.getArrayType(), "test3");
		assertEquals(data2.getEntryCount(), 2);

		ParameterNameValue param = null;
		Iterator<ParameterNameValue> it = data2.getAlgParams().iterator();
		param = it.next();
		assertEquals(param.getName(), "an1");
		assertEquals(param.getValueText(), "av1");

		it = data2.getSummaryParams().iterator();
		param = (ParameterNameValue)it.next();
		assertEquals(param.getName(), "sn1");
		assertEquals(param.getValueText(), "sv1");

		psQuantData = data2.getQuantificationEntry(0);
		assertEquals(psQuantData.getQuantification(), 10.0f, 0.0001f);
		assertEquals(psQuantData.getName(), "abc");
		psQuantData = data2.getQuantificationEntry(1);
		assertEquals(psQuantData.getQuantification(), 20.0f, 0.0001f);
		assertEquals(psQuantData.getName(), "xyz");
	}

	public void testWriteId() throws Exception {
		List<ParameterNameValue> params = new ArrayList<ParameterNameValue>();
		CHPQuantificationData data = new CHPQuantificationData("CHP_quantification_file_id");
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

		CHPQuantificationFileWriter writer = new CHPQuantificationFileWriter(data);
		ProbeSetQuantificationData psQuantData = new ProbeSetQuantificationData();

		writer.seekToDataSet();
		psQuantData.setId(10);
		psQuantData.setQuantification(10.0f);
		writer.writeEntry(psQuantData);
		psQuantData.setId(20);
		psQuantData.setQuantification(20.0f);
		writer.writeEntry(psQuantData);
		writer.close();
		assertTrue(true);

		CHPQuantificationData data2 = new CHPQuantificationData();
		CHPQuantificationFileReader reader = new CHPQuantificationFileReader();
		reader.setFilename("CHP_quantification_file_id");
		reader.read(data2);

		assertEquals(data2.getAlgName(), "sig");
		assertEquals(data2.getAlgVersion(), "1.0");
		assertEquals(data2.getArrayType(), "test3");
		assertEquals(data2.getEntryCount(), 2);

		ParameterNameValue param = null;
		Iterator it = data2.getAlgParams().iterator();
		param = (ParameterNameValue)it.next();
		assertEquals(param.getName(), "an1");
		assertEquals(param.getValueText(), "av1");

		it = data2.getSummaryParams().iterator();
		param = (ParameterNameValue)it.next();
		assertEquals(param.getName(), "sn1");
		assertEquals(param.getValueText(), "sv1");

		psQuantData = data2.getQuantificationEntry(0);
		assertEquals(psQuantData.getQuantification(), 10.0f, 0.0001f);
		assertEquals(psQuantData.getId(), 10);
		assertEquals(psQuantData.getName(), "");
		psQuantData = data2.getQuantificationEntry(1);
		assertEquals(psQuantData.getQuantification(), 20.0f, 0.0001f);
		assertEquals(psQuantData.getId(), 20);
		assertEquals(psQuantData.getName(), "");
	}
}
