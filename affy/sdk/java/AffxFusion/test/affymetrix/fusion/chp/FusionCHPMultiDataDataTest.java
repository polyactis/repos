/*
 * FusionCHPMultiDataDataTest.java
 *
 * Created on October 31, 2006, 10:37 AM
 *
 * To change this template, choose Tools | Options and locate the template under
 * the Source Creation and Management node. Right-click the template and choose
 * Open. You can then make changes to the template in the Source Editor.
 */

package affymetrix.fusion.chp;

import java.io.File;
import java.util.List;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import affymetrix.calvin.data.ProbeSetMultiDataCopyNumberData;
import affymetrix.calvin.data.ProbeSetMultiDataExpressionData;
import affymetrix.calvin.data.ProbeSetMultiDataGenotypeData;
import affymetrix.calvin.data.CHPMultiDataData.MultiDataType;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.calvin.parameter.ParameterNameValue.ParameterType;

/**
 * 
 * @author ljevon
 */
public class FusionCHPMultiDataDataTest extends TestCase {

	/** Creates a new instance of FusionCHPMultiDataDataTest */
	public FusionCHPMultiDataDataTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() throws Exception {
		TestSuite suite = new TestSuite(FusionCHPMultiDataDataTest.class);

		FusionCHPMultiDataData.registerReader();
		FusionCHPLegacyData.registerReader();

		return suite;
	}

	public void testReadNonMultiData() throws Exception {
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_expression_file");
		FusionCHPData chp = FusionCHPDataReg.read(f.getCanonicalPath());
		assertNotNull(chp);
		FusionCHPMultiDataData dataChp = FusionCHPMultiDataData.fromBase(chp);
		assertNull(dataChp);
	}

	public void testGenotype() throws Exception {
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_MultiData_file");
		FusionCHPData chp = FusionCHPDataReg.read(f.getCanonicalPath());
		assertNotNull(chp);
		FusionCHPMultiDataData dataChp = FusionCHPMultiDataData.fromBase(chp);
		assertNotNull(dataChp);

		assertEquals(dataChp.getAlgName(), "sig");
		assertEquals(dataChp.getAlgVersion(), "1.0");
		assertEquals(dataChp.getArrayType(), "test3");
		assertEquals(dataChp.getEntryCount(MultiDataType.ExpressionMultiDataType), 1);
		assertEquals(dataChp.getEntryCount(MultiDataType.GenotypeMultiDataType), 2);
		assertEquals(dataChp.getEntryCount(MultiDataType.ExpressionControlMultiDataType), 3);
		assertEquals(dataChp.getEntryCount(MultiDataType.GenotypeControlMultiDataType), 1);

		ProbeSetMultiDataExpressionData ex;
		ProbeSetMultiDataGenotypeData gn;
		ParameterNameValue param;

		// expression
		assertEquals(dataChp.getNumMetricColumns(MultiDataType.ExpressionMultiDataType), 0);
		ex = dataChp.getExpressionEntry(MultiDataType.ExpressionMultiDataType, 0);
		assertEquals(ex.getQuantification(), 10.0f, 0.0001f);
		assertEquals(ex.getName(), "ex1");

		// genotype
		assertEquals(dataChp.getNumMetricColumns(MultiDataType.GenotypeMultiDataType), 2);
		assertEquals(dataChp.getMetricColumnName(MultiDataType.GenotypeMultiDataType, 0), "int");
		assertEquals(dataChp.getMetricColumnName(MultiDataType.GenotypeMultiDataType, 1), "float");

		gn = dataChp.getGenotypeEntry(MultiDataType.GenotypeMultiDataType, 0);
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

		gn = dataChp.getGenotypeEntry(MultiDataType.GenotypeMultiDataType, 1);
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
		assertEquals(dataChp.getNumMetricColumns(MultiDataType.GenotypeControlMultiDataType), 0);
		gn = dataChp.getGenotypeEntry(MultiDataType.GenotypeControlMultiDataType, 0);
		assertEquals(gn.getCall(), 2);
		assertEquals(gn.getConfidence(), 22.0f, 0.0001f);
		assertEquals(gn.getMetrics().size(), 0);
		assertEquals(gn.getName(), "gc1");

		// expression control
		assertEquals(dataChp.getNumMetricColumns(MultiDataType.ExpressionControlMultiDataType), 0);
		ex = dataChp.getExpressionEntry(MultiDataType.ExpressionControlMultiDataType, 0);
		assertEquals(ex.getQuantification(), 20.0f, 0.0001f);
		assertEquals(ex.getName(), "ec1");
		ex = dataChp.getExpressionEntry(MultiDataType.ExpressionControlMultiDataType, 1);
		assertEquals(ex.getQuantification(), 30.0f, 0.0001f);
		assertEquals(ex.getName(), "ec2");
		ex = dataChp.getExpressionEntry(MultiDataType.ExpressionControlMultiDataType, 2);
		assertEquals(ex.getQuantification(), 40.0f, 0.0001f);
		assertEquals(ex.getName(), "ec3");
	}

	public void testCN() throws Exception {
		File f = new File("..\\..\\..\\..\\calvin_files\\parsers\\data\\CHP_MultiData_file_cn");
		FusionCHPData chp = FusionCHPDataReg.read(f.getCanonicalPath());
		assertNotNull(chp);
		FusionCHPMultiDataData cnChp = FusionCHPMultiDataData.fromBase(chp);
		assertNotNull(cnChp);

		ProbeSetMultiDataCopyNumberData e;
		ParameterNameValue param;

		assertEquals(cnChp.getAlgName(), "sig");
		assertEquals(cnChp.getAlgVersion(), "1.0");
		assertEquals(cnChp.getArrayType(), "test3");
		assertEquals(cnChp.getEntryCount(MultiDataType.CopyNumberMultiDataType), 2);
		assertEquals(cnChp.getEntryCount(MultiDataType.GenotypeMultiDataType), 0);
		assertEquals(cnChp.getEntryCount(MultiDataType.ExpressionMultiDataType), 0);
		assertEquals(cnChp.getEntryCount(MultiDataType.ExpressionControlMultiDataType), 0);
		assertEquals(cnChp.getEntryCount(MultiDataType.GenotypeControlMultiDataType), 0);

		List<ParameterNameValue> p = cnChp.getAlgParams();
		param = p.get(0);
		assertEquals(param.getName(), "an1");
		assertEquals(param.getValueText(), "av1");

		p = cnChp.getSummaryParams();
		param = p.get(0);
		assertEquals(param.getName(), "sn1");
		assertEquals(param.getValueText(), "sv1");

		e = cnChp.getCopyNumberEntry(MultiDataType.CopyNumberMultiDataType, 0);
		assertEquals(e.getChr(), (byte)10);
		assertEquals(e.getPosition(), 11);
		assertEquals(e.getName(), "abc");
		e = cnChp.getCopyNumberEntry(MultiDataType.CopyNumberMultiDataType, 1);
		assertEquals(e.getChr(), (byte)20);
		assertEquals(e.getPosition(), 21);
		assertEquals(e.getName(), "xyz");
	}

}
