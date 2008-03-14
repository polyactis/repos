package affymetrix.calvin.writers;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import junit.framework.TestCase;
import affymetrix.calvin.data.CHPMultiDataData;
import affymetrix.calvin.data.ColumnInfo;
import affymetrix.calvin.data.FloatColumn;
import affymetrix.calvin.data.IntColumn;
import affymetrix.calvin.data.ProbeSetMultiDataCopyNumberData;
import affymetrix.calvin.data.ProbeSetMultiDataCytoRegionData;
import affymetrix.calvin.data.ProbeSetMultiDataExpressionData;
import affymetrix.calvin.data.ProbeSetMultiDataGenotypeData;
import affymetrix.calvin.data.CHPMultiDataData.MultiDataType;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.parameter.ParameterNameValue.ParameterType;
import affymetrix.calvin.parsers.CHPMultiDataFileReader;
import affymetrix.portability.UInt;

public class CHPMultiDataFileWriterTest extends TestCase {

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	protected void setUp() throws Exception {
		super.setUp();
	}

	public void testWriteGeno() throws Exception {
		String filename = "CHP_MultiData_geno";
		CHPMultiDataData data = new CHPMultiDataData(filename);
		data.setAlgName("sig");
		data.setAlgVersion("1.0");
		data.setArrayType("test3");
		data.setEntryCount(MultiDataType.GenotypeMultiDataType, 2, 10);

		List<ParameterNameValue> params1 = new ArrayList<ParameterNameValue>();
		ParameterNameValue param1 = new ParameterNameValue();
		param1.setName("an1");
		param1.setValueText("av1");
		params1.add(param1);
		data.addAlgParams(params1);

		List<ParameterNameValue> params2 = new ArrayList<ParameterNameValue>();
		ParameterNameValue param2 = new ParameterNameValue();
		param2.setName("sn1");
		param2.setValueText("sv1");
		params2.add(param2);
		data.addSummaryParams(params2);

		CHPMultiDataFileWriter writer = new CHPMultiDataFileWriter(data);
		ProbeSetMultiDataGenotypeData e = new ProbeSetMultiDataGenotypeData();
		writer.seekToDataSet(MultiDataType.GenotypeMultiDataType);
		e.setName("abc");
		e.setConfidence(10.0f);
		writer.writeEntry(e);
		e.setName("xyz");
		e.setConfidence(20.0f);
		writer.writeEntry(e);
		writer.close();
		assertTrue(true);

		CHPMultiDataData data2 = new CHPMultiDataData();
		CHPMultiDataFileReader reader = new CHPMultiDataFileReader();
		reader.setFilename(filename);
		reader.read(data2);

		assertEquals(data2.getAlgName(), "sig");
		assertEquals(data2.getAlgVersion(), "1.0");
		assertEquals(data2.getArrayType(), "test3");
		assertEquals(data2.getEntryCount(MultiDataType.GenotypeMultiDataType), 2);
		assertEquals(data2.getEntryCount(MultiDataType.ExpressionMultiDataType), 0);
		assertEquals(data2.getEntryCount(MultiDataType.ExpressionControlMultiDataType), 0);
		assertEquals(data2.getEntryCount(MultiDataType.GenotypeControlMultiDataType), 0);

		List<ParameterNameValue> p = data2.getAlgParams();
		Iterator<ParameterNameValue> it = p.iterator();
		param1 = it.next();
		assertEquals(param1.getName(), "an1");
		assertEquals(param1.getValueText(), "av1");

		p = data2.getSummaryParams();
		it = p.iterator();
		param2 = it.next();
		assertEquals(param2.getName(), "sn1");
		assertEquals(param2.getValueText(), "sv1");

		e = data2.getGenotypeEntry(MultiDataType.GenotypeMultiDataType, 0);
		assertEquals(e.getConfidence(), 10.0f, 0.0001f);
		assertEquals(e.getName(), "abc");
		e = data2.getGenotypeEntry(MultiDataType.GenotypeMultiDataType, 1);
		assertEquals(e.getConfidence(), 20.0f, 0.0001f);
		assertEquals(e.getName(), "xyz");

	}

	public void testWriteCN() throws Exception {
		String filename = "CHP_MultiData_cn";
		CHPMultiDataData data = new CHPMultiDataData(filename);

		data.setAlgName("sig");
		data.setAlgVersion("1.0");
		data.setArrayType("test3");
		data.setEntryCount(MultiDataType.CopyNumberMultiDataType, 2, 10);
		data.setEntryCount(MultiDataType.CytoMultiDataType, 2, 10);

		List<ParameterNameValue> params1 = new ArrayList<ParameterNameValue>();
		ParameterNameValue param1 = new ParameterNameValue();
		param1.setName("an1");
		param1.setValueText("av1");
		params1.add(param1);
		data.addAlgParams(params1);

		List<ParameterNameValue> params2 = new ArrayList<ParameterNameValue>();
		ParameterNameValue param2 = new ParameterNameValue();
		params2.clear();
		param2.setName("sn1");
		param2.setValueText("sv1");
		params2.add(param2);
		data.addSummaryParams(params2);

		CHPMultiDataFileWriter writer = new CHPMultiDataFileWriter(data);
		ProbeSetMultiDataCopyNumberData e = new ProbeSetMultiDataCopyNumberData();
		ProbeSetMultiDataCytoRegionData c = new ProbeSetMultiDataCytoRegionData();

		writer.seekToDataSet(MultiDataType.CopyNumberMultiDataType);
		e.setName("abc");
		e.setChr((byte)10);
		e.setPosition(11);
		writer.writeEntry(e);
		e.setName("xyz");
		e.setChr((byte)20);
		e.setPosition(21);
		writer.writeEntry(e);

		writer.seekToDataSet(MultiDataType.CytoMultiDataType);
		c.setName("abc");
		c.setCall((byte)10);
		c.setConfidence(11f);
		c.setChr((byte)10);
		c.setStartPosition(new UInt(10));
		c.setStopPosition(new UInt(11));
		writer.writeEntry(c);

		c.setName("xyz");
		c.setCall((byte)20);
		c.setConfidence(21f);
		c.setChr((byte)20);
		c.setStartPosition(new UInt(20));
		c.setStopPosition(new UInt(21));
		writer.writeEntry(c);
		writer.close();
		assertTrue(true);

		CHPMultiDataData data2 = new CHPMultiDataData();
		CHPMultiDataFileReader reader = new CHPMultiDataFileReader();
		reader.setFilename(filename);
		reader.read(data2);

		assertEquals(data2.getAlgName(), "sig");
		assertEquals(data2.getAlgVersion(), "1.0");
		assertEquals(data2.getArrayType(), "test3");
		assertEquals(data2.getEntryCount(MultiDataType.CopyNumberMultiDataType), 2);
		assertEquals(data2.getEntryCount(MultiDataType.GenotypeMultiDataType), 0);
		assertEquals(data2.getEntryCount(MultiDataType.ExpressionMultiDataType), 0);
		assertEquals(data2.getEntryCount(MultiDataType.ExpressionControlMultiDataType), 0);
		assertEquals(data2.getEntryCount(MultiDataType.GenotypeControlMultiDataType), 0);

		List<ParameterNameValue> p = data2.getAlgParams();
		Iterator<ParameterNameValue> it = p.iterator();
		param1 = it.next();
		assertEquals(param1.getName(), "an1");
		assertEquals(param1.getValueText(), "av1");

		p = data2.getSummaryParams();
		it = p.iterator();
		param2 = it.next();
		assertEquals(param2.getName(), "sn1");
		assertEquals(param2.getValueText(), "sv1");

		e = data2.getCopyNumberEntry(MultiDataType.CopyNumberMultiDataType, 0);
		assertEquals(e.getChr(), 10);
		assertEquals(e.getPosition(), 11);
		assertEquals(e.getName(), "abc");
		e = data2.getCopyNumberEntry(MultiDataType.CopyNumberMultiDataType, 1);
		assertEquals(e.getChr(), 20);
		assertEquals(e.getPosition(), 21);
		assertEquals(e.getName(), "xyz");

		c = data2.getCytoEntry(MultiDataType.CytoMultiDataType, 0);
		assertEquals(c.getCall(), 10);
		assertEquals(c.getChr(), 10);
		assertEquals(c.getStartPosition().toLong(), 10);
		assertEquals(c.getStopPosition().toLong(), 11);
		assertEquals(c.getConfidence(), 11.0f, 0.0001f);
		assertEquals(c.getName(), "abc");

		c = data2.getCytoEntry(MultiDataType.CytoMultiDataType, 1);
		assertEquals(c.getCall(), 20);
		assertEquals(c.getChr(), 20);
		assertEquals(c.getStartPosition().toLong(), 20);
		assertEquals(c.getStopPosition().toLong(), 21);
		assertEquals(c.getConfidence(), 21.0f, 0.0001f);
		assertEquals(c.getName(), "xyz");
	}

	public void testWriteExp() throws Exception {
		String filename = "CHP_MultiData_exp";
		CHPMultiDataData data = new CHPMultiDataData(filename);
		data.setAlgName("sig");
		data.setAlgVersion("1.0");
		data.setArrayType("test3");
		data.setEntryCount(MultiDataType.ExpressionMultiDataType, 2, 10);

		CHPMultiDataFileWriter writer = new CHPMultiDataFileWriter(data);
		ProbeSetMultiDataExpressionData e = new ProbeSetMultiDataExpressionData();
		writer.seekToDataSet(MultiDataType.ExpressionMultiDataType);
		e.setName("abc");
		e.setQuantification(10f);
		writer.writeEntry(e);
		e.setName("xyz");
		e.setQuantification(20f);
		writer.writeEntry(e);
		writer.close();
		assertTrue(true);

		CHPMultiDataData data2 = new CHPMultiDataData();
		CHPMultiDataFileReader reader = new CHPMultiDataFileReader();
		reader.setFilename(filename);
		reader.read(data2);
		assertEquals(data2.getAlgName(), "sig");
		assertEquals(data2.getAlgVersion(), "1.0");
		assertEquals(data2.getArrayType(), "test3");
		assertEquals(data2.getEntryCount(MultiDataType.ExpressionMultiDataType), 2);
		assertEquals(data2.getEntryCount(MultiDataType.GenotypeMultiDataType), 0);
		assertEquals(data2.getEntryCount(MultiDataType.ExpressionControlMultiDataType), 0);
		assertEquals(data2.getEntryCount(MultiDataType.GenotypeControlMultiDataType), 0);
		e = data2.getExpressionEntry(MultiDataType.ExpressionMultiDataType, 0);
		assertEquals(e.getQuantification(), 10.0f, 0.0001f);
		assertEquals(e.getName(), "abc");
		e = data2.getExpressionEntry(MultiDataType.ExpressionMultiDataType, 1);
		assertEquals(e.getQuantification(), 20.0f, 0.0001f);
		assertEquals(e.getName(), "xyz");
	}

	public void testWriteAll() throws Exception {
		String filename = "CHP_MultiData_all";
		CHPMultiDataData data = new CHPMultiDataData(filename);

		List<ColumnInfo> cols = new ArrayList<ColumnInfo>();
		IntColumn icol = new IntColumn("int");
		cols.add(icol);
		FloatColumn fcol = new FloatColumn("float");
		cols.add(fcol);

		data.setAlgName("sig");
		data.setAlgVersion("1.0");
		data.setArrayType("test3");
		data.setEntryCount(MultiDataType.ExpressionMultiDataType, 1, 10);
		data.setEntryCount(MultiDataType.GenotypeMultiDataType, 2, 10, cols);
		data.setEntryCount(MultiDataType.ExpressionControlMultiDataType, 3, 10);
		data.setEntryCount(MultiDataType.GenotypeControlMultiDataType, 1, 10);

		CHPMultiDataFileWriter writer = new CHPMultiDataFileWriter(data);
		ProbeSetMultiDataExpressionData ex = new ProbeSetMultiDataExpressionData();
		ProbeSetMultiDataGenotypeData gn = new ProbeSetMultiDataGenotypeData();

		// expression
		writer.seekToDataSet(MultiDataType.ExpressionMultiDataType);
		ex.setName("ex1");
		ex.setQuantification(10f);
		writer.writeEntry(ex);

		// expression control
		writer.seekToDataSet(MultiDataType.ExpressionControlMultiDataType);
		ex.setName("ec1");
		ex.setQuantification(20f);
		writer.writeEntry(ex);
		ex.setName("ec2");
		ex.setQuantification(30f);
		writer.writeEntry(ex);

		// genotype
		writer.seekToDataSet(MultiDataType.GenotypeMultiDataType);
		gn.setName("gn1");
		gn.setCall((byte)1);
		gn.setConfidence(11f);
		ParameterNameValue param1 = new ParameterNameValue();
		param1.setValueInt32(1);
		gn.addMetric(param1);
		ParameterNameValue param2 = new ParameterNameValue();
		param2.setValueFloat(2.0f);
		gn.addMetric(param2);
		writer.writeEntry(gn);

		gn.setName("gn2");
		gn.setCall((byte)2);
		gn.setConfidence(22f);
		gn.clearMetrics();
		ParameterNameValue param3 = new ParameterNameValue();
		param3.setValueInt32(2);
		gn.addMetric(param3);
		ParameterNameValue param4 = new ParameterNameValue();
		param4.setValueFloat(3.0f);
		gn.addMetric(param4);
		writer.writeEntry(gn);

		// genotype control
		writer.seekToDataSet(MultiDataType.GenotypeControlMultiDataType);
		gn.setName("gc1");
		gn.setCall((byte)2);
		gn.setConfidence(22f);
		gn.clearMetrics();
		writer.writeEntry(gn);

		// expression control
		writer.seekToDataSet(MultiDataType.ExpressionControlMultiDataType);
		ex.setName("ec3");
		ex.setQuantification(40f);
		writer.writeEntry(ex);
		writer.close();
		assertTrue(true);

		CHPMultiDataData data2 = new CHPMultiDataData();
		CHPMultiDataFileReader reader = new CHPMultiDataFileReader();
		reader.setFilename(filename);
		reader.read(data2);

		assertEquals(data2.getAlgName(), "sig");
		assertEquals(data2.getAlgVersion(), "1.0");
		assertEquals(data2.getArrayType(), "test3");
		assertEquals(data2.getEntryCount(MultiDataType.ExpressionMultiDataType), 1);
		assertEquals(data2.getEntryCount(MultiDataType.GenotypeMultiDataType), 2);
		assertEquals(data2.getEntryCount(MultiDataType.ExpressionControlMultiDataType), 3);
		assertEquals(data2.getEntryCount(MultiDataType.GenotypeControlMultiDataType), 1);

		// expression
		assertEquals(data2.getNumMetricColumns(MultiDataType.ExpressionMultiDataType), 0);
		ex = data2.getExpressionEntry(MultiDataType.ExpressionMultiDataType, 0);
		assertEquals(ex.getQuantification(), 10.0f, 0.0001f);
		assertEquals(ex.getName(), "ex1");

		// genotype
		assertEquals(data2.getNumMetricColumns(MultiDataType.GenotypeMultiDataType), 2);
		assertEquals(data2.getMetricColumnName(MultiDataType.GenotypeMultiDataType, 0), "int");
		assertEquals(data2.getMetricColumnName(MultiDataType.GenotypeMultiDataType, 1), "float");

		ParameterNameValue param = null;
		gn = data2.getGenotypeEntry(MultiDataType.GenotypeMultiDataType, 0);
		assertEquals(gn.getCall(), 1);
		assertEquals(gn.getConfidence(), 11.0f, 0.0001f);
		assertEquals(gn.getName(), "gn1");
		assertEquals(gn.getMetrics().size(), 2);
		param = gn.getMetrics().get(0);
		assertEquals(param.getParameterType(), ParameterType.Int32Type);
		assertEquals(param.getValueInt32(), 1);
		param = gn.getMetrics().get(1);
		assertEquals(param.getParameterType(), ParameterType.FloatType);
		assertEquals(param.getValueFloat(), 2.0f, 0.00001f);

		gn = data2.getGenotypeEntry(MultiDataType.GenotypeMultiDataType, 1);
		assertEquals(gn.getCall(), 2);
		assertEquals(gn.getConfidence(), 22.0f, 0.0001f);
		assertEquals(gn.getName(), "gn2");
		assertEquals(gn.getMetrics().size(), 2);
		param = gn.getMetrics().get(0);
		assertEquals(param.getParameterType(), ParameterType.Int32Type);
		assertEquals(param.getValueInt32(), 2);
		param = gn.getMetrics().get(1);
		assertEquals(param.getParameterType(), ParameterType.FloatType);
		assertEquals(param.getValueFloat(), 3.0f, 0.00001f);

		// genotype control
		assertEquals(data2.getNumMetricColumns(MultiDataType.GenotypeControlMultiDataType), 0);
		gn = data2.getGenotypeEntry(MultiDataType.GenotypeControlMultiDataType, 0);
		assertEquals(gn.getCall(), 2);
		assertEquals(gn.getConfidence(), 22.0f, 0.0001f);
		assertEquals(gn.getMetrics().size(), 0);
		assertEquals(gn.getName(), "gc1");

		// expression control
		assertEquals(data2.getNumMetricColumns(MultiDataType.ExpressionControlMultiDataType), 0);
		ex = data2.getExpressionEntry(MultiDataType.ExpressionControlMultiDataType, 0);
		assertEquals(ex.getQuantification(), 20.0f, 0.0001f);
		assertEquals(ex.getName(), "ec1");
		ex = data2.getExpressionEntry(MultiDataType.ExpressionControlMultiDataType, 1);
		assertEquals(ex.getQuantification(), 30.0f, 0.0001f);
		assertEquals(ex.getName(), "ec2");
		ex = data2.getExpressionEntry(MultiDataType.ExpressionControlMultiDataType, 2);
		assertEquals(ex.getQuantification(), 40.0f, 0.0001f);
		assertEquals(ex.getName(), "ec3");
	}
}
