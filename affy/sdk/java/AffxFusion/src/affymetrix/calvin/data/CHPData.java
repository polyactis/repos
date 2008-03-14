/////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 2.1 of the License,
// or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
/////////////////////////////////////////////////////////////////

package affymetrix.calvin.data;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.data.ColumnInfo.DataSetColumnTypes;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.AffymetrixParameterConsts;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.portability.UByte;
import affymetrix.portability.UShort;

/** Defines a base class for older CHP file data. */
public class CHPData extends ChpDataBase {

	/** The id for the expression CHP files. */
	public static final String CHP_EXPRESSION_ASSAY_TYPE = "affymetrix-expression-probeset-analysis";

	/** The id for the expression CHP file data group. */
	public static final String CHP_EXPR_GROUP = "Expression Results";

	/** The id for the resequencing CHP files. */
	public static final String CHP_RESEQUENCING_ASSAY_TYPE = "affymetrix-resequencing-probeset-analysis";

	/** The id for the resequencing CHP file data group. */
	public static final String CHP_RESEQ_GROUP = "Resequencing Results";

	/** The id for the genotyping CHP files. */
	public static final String CHP_GENOTYPING_ASSAY_TYPE = "affymetrix-genotyping-probeset-analysis";

	/** The id for the genotyping CHP file data group. */
	public static final String CHP_GENO_GROUP = "Genotyping Results";

	/** The id for the universal CHP files. */
	public static final String CHP_UNIVERSAL_ASSAY_TYPE = "affymetrix-universal-probeset-analysis";

	/** The id for the universal CHP file data group. */
	public static final String CHP_UNIV_GROUP = "Universal Results";

	/** The id for the number of rows of features. */
	public static final String CHP_ROWS = "affymetrix-cel-rows";

	/** The id for the number of columns of features. */
	public static final String CHP_COLS = "affymetrix-cel-cols";

	/** The id for the prog ID. */
	public static final String CHP_PROGID = "affymetrix-progid";

	/** The id for the parent cel file. */
	public static final String CHP_PARENT_CELL = "affymetrix-parent-celfile";

	/** The group name for the background zone group. */
	public static final String CHP_BG_ZONE_GROUP = "Background Zone Data";

	/** The group name for the force call group (for resequencing only). */
	public static final String CHP_RESEQ_FORCE_CALL_GROUP = "Force Call Data";

	/** The group name for the orig call group (for resequencing only). */
	public static final String CHP_RESEQ_ORIG_CALL_GROUP = "Orig Call Data";

	// Constant column names.
	private static final String CallColName = "Call";

	private static final String ScoreColName = "Score";

	private static final String BackgroundColName = "Background";

	private static final String ConfidenceColName = "Confidence";

	private static final String RAS1ColName = "RAS1";

	private static final String RAS2ColName = "RAS2";

	private static final String AAColName = "AA Call p-value";

	private static final String ABColName = "AB Call p-value";

	private static final String BBColName = "BB Call p-value";

	private static final String NoCallColName = "No Call p-value";

	private static final String ProbeSetNameColName = "Probe Set Name";

	private static final String DetectionColName = "Detection";

	private static final String DetectionPValueColName = "Detection p-value";

	private static final String SignalColName = "Signal";

	private static final String NumberPairsColName = "Number of Pairs";

	private static final String NumberPairsUsedColName = "Number of Pairs Used";

	private static final String ChangeColName = "Change";

	private static final String ChangePValueColName = "Change p-value";

	private static final String SignalLogRatioColName = "Signal Log Ratio";

	private static final String SignalLogRatioLowColName = "Signal Log Ratio Low";

	private static final String SignalLogRatioHighColName = "Signal Log Ratio High";

	private static final String CommonPairsColName = "Common Pairs";

	private static final String CenterXColName = "Center X";

	private static final String CenterYColName = "Center Y";

	private static final String SmoothFactorColName = "Smooth Factor";

	private static final String PositionColName = "Position";

	private static final String ReasonColName = "Reason";

	private static final String ForceCallColName = "Force Call";

	private static final String OriginalCallColName = "Original Call";

	/* ! Flag indicating if the probe set names were stored in wide character format. */
	private boolean wideProbeSetNames = false;

	/** keep rows from being read from the header all the time */
	private int cachedRows = -1;

	/** keep cols from being read from the header all the time */
	private int cachedCols = -1;

	/** expression entries DataSet */
	private DataSet entriesExp = null;

	/** genotyping entries DataSet */
	private DataSet entriesGeno = null;

	/** universal entries DataSet */
	private DataSet entriesUniv = null;

	/** resequencing entries DataSet */
	private DataSet entriesReseq = null;

	/** chp background zones DataSet */
	private DataSet bgZones = null;

	/** chp force call DataSet */
	private DataSet forceSet = null;

	/** chp orig DataSet */
	private DataSet origSet = null;

	/* ! The maximum length of a probe set name. */
	int maxProbeSetName = -1;

	/** Constructor */
	public CHPData() {
	}

	/**
	 * Constructor with file and type.
	 * 
	 * @param filename
	 *          The name of the CHP file.
	 * @param assayType
	 *          The type of data in the CHP file.
	 */
	public CHPData(String filename, String assayType) {

		this();
		setFilename(filename);

		String groupName;
		if (assayType.equals(CHP_EXPRESSION_ASSAY_TYPE)) {
			groupName = CHP_EXPR_GROUP;
		}
		else if (assayType.equals(CHP_RESEQUENCING_ASSAY_TYPE)) {
			groupName = CHP_RESEQ_GROUP;
		}
		else if (assayType.equals(CHP_GENOTYPING_ASSAY_TYPE)) {
			groupName = CHP_GENO_GROUP;
		}
		else if (assayType.equals(CHP_UNIVERSAL_ASSAY_TYPE)) {
			groupName = CHP_UNIV_GROUP;
		}
		else {
			return;
		}

		DataGroupHeader dcHdr = new DataGroupHeader(groupName);
		genericData.getHeader().addDataGroupHdr(dcHdr);
		genericData.getHeader().getGenericDataHdr().setFileTypeId(assayType);
		DataGroupHeader dcHdrBg = new DataGroupHeader(CHP_BG_ZONE_GROUP);
		genericData.getHeader().addDataGroupHdr(dcHdrBg);

		// Now add the force and orig sets for reseq designs.
		if (assayType.equals(CHP_RESEQUENCING_ASSAY_TYPE)) {
			DataGroupHeader dcHdrForce = new DataGroupHeader(CHP_RESEQ_FORCE_CALL_GROUP);
			genericData.getHeader().addDataGroupHdr(dcHdrForce);
			DataGroupHeader dcHdrOrig = new DataGroupHeader(CHP_RESEQ_ORIG_CALL_GROUP);
			genericData.getHeader().addDataGroupHdr(dcHdrOrig);
		}
	}

	void addColumns(DataSetHeader hdr, boolean hasCompData) {
		if (genericData.getHeader().getGenericDataHdr().getFileTypeId() == CHP_EXPRESSION_ASSAY_TYPE) {
			addExprColumns(hdr, hasCompData);
		}

		if (genericData.getHeader().getGenericDataHdr().getFileTypeId() == CHP_GENOTYPING_ASSAY_TYPE) {
			addGenoColumns(hdr);
		}
		else if (genericData.getHeader().getGenericDataHdr().getFileTypeId() == CHP_UNIVERSAL_ASSAY_TYPE) {
			addUnivColumns(hdr);
		}
		else if (genericData.getHeader().getGenericDataHdr().getFileTypeId() == CHP_RESEQUENCING_ASSAY_TYPE) {
			addReseqColumns(hdr);
		}
	}

	protected void addExprColumns(DataSetHeader hdr, boolean hasCompData) {
		hdr.setName(CHP_EXPR_GROUP);
		// Probeset name - string
		hdr.addAsciiColumn(ProbeSetNameColName, maxProbeSetName);
		// Detection - unsigned char
		hdr.addUByteColumn(DetectionColName);
		// Detection p-value - float
		hdr.addFloatColumn(DetectionPValueColName);
		// Signal - float
		hdr.addFloatColumn(SignalColName);
		// Number of pairs - unsigned short
		hdr.addUShortColumn(NumberPairsColName);
		// Number of pairs used - unsigned short
		hdr.addUShortColumn(NumberPairsUsedColName);

		if (hasCompData == true) {
			// Change - unsigned char
			hdr.addUByteColumn(ChangeColName);
			// Change p-value - float
			hdr.addFloatColumn(ChangePValueColName);
			// Signal Log Ratio - float
			hdr.addFloatColumn(SignalLogRatioColName);
			// Signal Log Ratio Low - float
			hdr.addFloatColumn(SignalLogRatioLowColName);
			// Signal Log Ratio High - float
			hdr.addFloatColumn(SignalLogRatioHighColName);
			// Common Pairs - unsigned short
			hdr.addUShortColumn(CommonPairsColName);
		}
	}

	void addGenoColumns(DataSetHeader hdr) {
		hdr.setName(CHP_GENO_GROUP);
		// Probeset name - string
		hdr.addAsciiColumn(ProbeSetNameColName, maxProbeSetName);
		// Call - unsigned char
		hdr.addUByteColumn(CallColName);
		// Confidence - float
		hdr.addFloatColumn(ConfidenceColName);
		// RAS1 - float
		hdr.addFloatColumn(RAS1ColName);
		// RAS2 - float
		hdr.addFloatColumn(RAS2ColName);
		// AA Call - float
		hdr.addFloatColumn(AAColName);
		// AB Call - float
		hdr.addFloatColumn(ABColName);
		// BB Call - float
		hdr.addFloatColumn(BBColName);
		// No Call - float
		hdr.addFloatColumn(NoCallColName);
	}

	void addReseqColumns(DataSetHeader hdr) {
		hdr.setName(CHP_RESEQ_GROUP);
		// call - char
		hdr.addByteColumn(CallColName);
		// Score - float
		hdr.addFloatColumn(ScoreColName);
	}

	void addUnivColumns(DataSetHeader hdr) {
		hdr.setName(CHP_UNIV_GROUP);
		// Background - float
		hdr.addFloatColumn(BackgroundColName);
	}

	/** Clears the members. */
	public void clear() {
		entriesGeno = null;
		entriesExp = null;
		entriesUniv = null;
		entriesReseq = null;
		bgZones = null;
		forceSet = null;
		origSet = null;
		genericData.getHeader().clear();
		cachedRows = -1;
		cachedCols = -1;
	}

	/** Gets the algorithm name */
	public String getAlgName() {
		return getStringFromGenericHdr(AffymetrixParameterConsts.ALGORITHM_NAME_PARAM_NAME);
	}

	/** Gets a single algorithm parameter by name. */
	public ParameterNameValue getAlgParam(String tag) {
		GenericDataHeader hdr = genericData.getHeader().getGenericDataHdr();
		ParameterNameValue type = hdr.findNameValParam(AffymetrixParameterConsts.ALGORITHM_PARAM_NAME_PREFIX + tag);
		if (type != null) {
			type.setName(tag);
		}
		return type;
	}

	/** Gets the alg parameters */
	public List<ParameterNameValue> getAlgParams() {
		return getParams(AffymetrixParameterConsts.ALGORITHM_PARAM_NAME_PREFIX);
	}

	/** Gets the algorithm version. */
	public String getAlgVersion() {
		return getStringFromGenericHdr(AffymetrixParameterConsts.ALG_VERSION_PARAM_NAME);
	}

	/** Gets the array type */
	public String getArrayType() {
		return getStringFromGenericHdr(AffymetrixParameterConsts.ARRAY_TYPE_PARAM_NAME);
	}

	/** Gets the assay type */
	public String getAssayType() {
		return genericData.getHeader().getGenericDataHdr().getFileTypeId();
	}

	/**
	 * Gets CHP background zone value
	 * 
	 * @param row
	 *          The row from which to start copying
	 * @param zone
	 *          The data object to be filled
	 */
	public void getBackgroundZone(int row, CHPBackgroundZone zone) throws UnsignedOutOfLimitsException {
		prepareBackgroundZoneDataSet();
		if ((bgZones != null) && bgZones.isOpen()) {
			float centerX = bgZones.getDataFloat(row, 0);
			zone.setCenterX(centerX);

			float centerY = bgZones.getDataFloat(row, 1);
			zone.setCenterY(centerY);

			float background = bgZones.getDataFloat(row, 2);
			zone.setBackground(background);

			float smoothFactor = bgZones.getDataFloat(row, 3);
			zone.setSmoothFactor(smoothFactor);
		}
	}

	/** Gets the number of background zones. */
	public int getBackgroundZoneCnt() {
		DataGroupHeader dcHdr = genericData.getHeader().getDataGroup(1);
		DataSetHeader dpHdr = dcHdr.getDataSet(0);
		return dpHdr.getRowCnt();
	}

	/** Gets the background zones. */
	public void getBackgroundZones(int row, int rowCnt, List<CHPBackgroundZone> zones)
			throws UnsignedOutOfLimitsException {
		zones.clear();
		for (int i = row; i < rowCnt; i++) {
			CHPBackgroundZone z = new CHPBackgroundZone();
			getBackgroundZone(i, z);
			zones.add(z);
		}
	}

	/** Gets a chip summary parameter by name */
	public ParameterNameValue getChipSum(String tag) {
		GenericDataHeader hdr = genericData.getHeader().getGenericDataHdr();
		ParameterNameValue type = hdr.findNameValParam(AffymetrixParameterConsts.CHIP_SUMMARY_PARAM_NAME_PREFIX + tag);
		if (type != null) {
			type.setName(tag);
		}
		return type;
	}

	/** Gets all the chip summary parameters */
	public List<ParameterNameValue> getChipSums() {
		List<ParameterNameValue> sumParams = new ArrayList<ParameterNameValue>();
		List<ParameterNameValue> allParams = genericData.getHeader().getGenericDataHdr().getNameValParams();
		for (int i = 0; i < allParams.size(); i++) {
			ParameterNameValue param = allParams.get(i);
			String name = param.getName();
			if (name.startsWith(AffymetrixParameterConsts.CHIP_SUMMARY_PARAM_NAME_PREFIX) == true) {
				ParameterNameValue sumParam = new ParameterNameValue(param);
				sumParam.setName(name.substring(AffymetrixParameterConsts.CHIP_SUMMARY_PARAM_NAME_PREFIX.length(), name
						.length()));
				sumParams.add(sumParam);
			}
		}
		return sumParams;
	}

	/** Gets the number of columns of features on the array. */
	public int getCols() {
		if (cachedCols == -1) {
			cachedCols = getInt32FromGenericHdr(CHP_COLS);
		}
		return cachedCols;
	}

	/** Gets the expression entry (probe set). */
	public void getEntry(int row, CHPExpressionEntry e) throws IOException, UnsignedOutOfLimitsException {
		prepareExprEntryDataSet();
		if ((entriesExp != null) && entriesExp.isOpen()) {
			int colIndex = 0;
			String probeSetName;
			if (wideProbeSetNames == false) {
				probeSetName = entriesExp.getDataString8(row, colIndex);
			}
			else {
				probeSetName = entriesExp.getDataString16(row, colIndex);
			}
			++colIndex;
			e.setProbeSetName(probeSetName);

			UByte detection = entriesExp.getDataUByte(row, colIndex);
			++colIndex;
			e.setDetection(detection);

			float detectionPValue = entriesExp.getDataFloat(row, colIndex);
			++colIndex;
			e.setDetectionPValue(detectionPValue);

			float signal = entriesExp.getDataFloat(row, colIndex);
			++colIndex;
			e.setSignal(signal);

			UShort numPairs = entriesExp.getDataUShort(row, colIndex);
			++colIndex;
			e.setNumPairs(numPairs);

			UShort numPairsUsed = entriesExp.getDataUShort(row, colIndex);
			++colIndex;
			e.setNumPairsUsed(numPairsUsed);

			e.setHasComparisonData(false);

			if (entriesExp.getCols() > colIndex) {
				e.setHasComparisonData(true);

				UByte change = entriesExp.getDataUByte(row, colIndex);
				++colIndex;
				e.setChange(change);

				float changePValue = entriesExp.getDataFloat(row, colIndex);
				++colIndex;
				e.setChangePValue(changePValue);

				float sigLogRatio = entriesExp.getDataFloat(row, colIndex);
				++colIndex;
				e.setSigLogRatio(sigLogRatio);

				float sigLogRatioLo = entriesExp.getDataFloat(row, colIndex);
				++colIndex;
				e.setSigLogRatioLo(sigLogRatioLo);

				float sigLogRatioHi = entriesExp.getDataFloat(row, colIndex);
				++colIndex;
				e.setSigLogRatioHi(sigLogRatioHi);

				UShort commonPairs = entriesExp.getDataUShort(row, colIndex);
				++colIndex;
				e.setCommonPairs(commonPairs);
			}
		}
	}

	/** Gets the expression entry (probe set). */
	public void getEntry(int row, CHPGenotypeEntry e) throws IOException, UnsignedOutOfLimitsException {
		prepareGenoEntryDataSet();
		if ((entriesGeno != null) && (entriesGeno.isOpen() == true)) {
			String probeSetName;
			if (wideProbeSetNames == false) {
				probeSetName = entriesGeno.getDataString8(row, 0);
			}
			else {
				probeSetName = entriesGeno.getDataString16(row, 0);
			}
			e.setProbeSetName(probeSetName);

			byte call = entriesGeno.getDataByte(row, 1);
			e.setCall(call);

			float confidence = entriesGeno.getDataFloat(row, 2);
			e.setConfidence(confidence);

			float ras1 = entriesGeno.getDataFloat(row, 3);
			e.setRAS1(ras1);

			float ras2 = entriesGeno.getDataFloat(row, 4);
			e.setRAS2(ras2);

			float aaCall = entriesGeno.getDataFloat(row, 5);
			e.setAACall(aaCall);

			float abCall = entriesGeno.getDataFloat(row, 6);
			e.setABCall(abCall);

			float bbCall = entriesGeno.getDataFloat(row, 7);
			e.setBBCall(bbCall);

			float noCall = entriesGeno.getDataFloat(row, 8);
			e.setNoCall(noCall);
		}
	}

	/** Gets the resequence entry. */
	public void getEntry(int row, CHPReseqEntry e) throws UnsignedOutOfLimitsException {
		prepareReseqEntryDataSet();
		if ((entriesReseq != null) && (entriesReseq.isOpen() == true)) {
			e.setCall(entriesReseq.getDataByte(row, 0));
			e.setScore(entriesReseq.getDataFloat(row, 1));
		}
	}

	/** Gets the expression entry (probe set). */
	public void getEntry(int row, CHPUniversalEntry e) throws UnsignedOutOfLimitsException {
		prepareUnivEntryDataSet();
		if ((entriesUniv != null) && (entriesUniv.isOpen() == true)) {
			float background = entriesUniv.getDataFloat(row, 0);
			e.setBackground(background);
		}
	}

	/** Gets the number of entries (probe sets) */
	public int getEntryCount() {
		DataGroupHeader dcHdr = genericData.getHeader().getDataGroup(0);
		DataSetHeader dpHdr = dcHdr.getDataSet(0);
		return dpHdr.getRowCnt();
	}

	/** Gets the file header */
	public FileHeader getFileHeader() {
		return genericData.getHeader();
	}

	/** The name of the CHP file. */
	public String getFilename() {
		return genericData.getHeader().getFilename();
	}

	/**
	 * Gets the force call value
	 * 
	 * @param row
	 *          The row index
	 * @param force
	 *          The data object to be filled
	 */
	public void getForceCall(int row, CHPReseqForceCall force) throws UnsignedOutOfLimitsException {
		prepareForceDataSet();
		if ((forceSet != null) && (forceSet.isOpen() == true)) {
			force.setPosition(forceSet.getDataInt(row, 0));
			force.setCall(forceSet.getDataByte(row, 1));
			force.setReason(forceSet.getDataByte(row, 2));
		}
	}

	/** Gets the number of force calls. */
	public int getForceCnt() {
		DataGroupHeader dcHdr = genericData.getHeader().getDataGroup(2);
		DataSetHeader dpHdr = dcHdr.getDataSet(0);
		return dpHdr.getRowCnt();
	}

	public void setForceCnt(int ln) {
		DataSetHeader dpHdr = new DataSetHeader();
		dpHdr.setRowCnt(ln);
		dpHdr.setName(CHP_RESEQ_FORCE_CALL_GROUP);

		// position - int
		dpHdr.addIntColumn(PositionColName);
		// call - byte
		dpHdr.addByteColumn(ForceCallColName);
		// reason - byte
		dpHdr.addByteColumn(ReasonColName);

		DataGroupHeader dcHdr = genericData.getHeader().getDataGroup(2);
		dcHdr.addDataSetHdr(dpHdr);
	}

	public void setOrigCnt(int ln) {
		DataSetHeader dpHdr = new DataSetHeader();
		dpHdr.setRowCnt(ln);
		dpHdr.setName(CHP_RESEQ_ORIG_CALL_GROUP);

		// position - int
		dpHdr.addIntColumn(PositionColName);
		// call - byte
		dpHdr.addByteColumn(OriginalCallColName);

		DataGroupHeader dcHdr = genericData.getHeader().getDataGroup(3);
		dcHdr.addDataSetHdr(dpHdr);
	}

	/** Gets the file data object. */
	public GenericData getGenericData() {
		return genericData;
	}

	/** Gets an integer parameter from the header. */
	private int getInt32FromGenericHdr(String name) {
		GenericDataHeader hdr = genericData.getHeader().getGenericDataHdr();
		ParameterNameValue paramType = hdr.findNameValParam(name);
		if (paramType != null) {
			return paramType.getValueInt32();
		}
		return 0;
	}

	/** Gets the files magic number. */
	public byte getMagic() {
		return genericData.getHeader().getMagicNumber();
	}

	/* ! The maximum length of a probe set name. */
	public int getMaxProbeSetName() {
		return maxProbeSetName;
	}

	/**
	 * Gets the original call value from the orig set.
	 * 
	 * @param row
	 *          The row index
	 * @param orig
	 *          The orginal call value.
	 */
	public void getOrigCall(int row, CHPReseqOrigCall orig) throws UnsignedOutOfLimitsException {
		prepareOrigDataSet();
		if ((origSet != null) && (origSet.isOpen() == true)) {
			orig.setPosition(origSet.getDataInt(row, 0));
			orig.setCall(origSet.getDataByte(row, 1));
		}
	}

	/** Gets the number of orig calls. */
	public int getOrigCnt() {
		DataGroupHeader dcHdr = genericData.getHeader().getDataGroup(3);
		DataSetHeader dpHdr = dcHdr.getDataSet(0);
		return dpHdr.getRowCnt();
	}

	/** Gets the name of the parent CEL file. */
	public String getParentCell() {
		return getStringFromGenericHdr(CHP_PARENT_CELL);
	}

	/** Gets the CHP file prog Id. */
	public String getProgId() {
		return getStringFromGenericHdr(CHP_PROGID);
	}

	/** Gets the number of rows of features on the array. */
	public int getRows() {
		if (cachedRows == -1) {
			cachedRows = getInt32FromGenericHdr(CHP_ROWS);
		}
		return cachedRows;
	}

	/** Gets the version in the file. */
	public byte getVersion() {
		return genericData.getHeader().getVersion();
	}

	/** Prepares the data set for the bg set. */
	private void prepareBackgroundZoneDataSet() {
		if (bgZones == null) {
			try {
				bgZones = genericData.getDataSet(CHP_BG_ZONE_GROUP, CHP_BG_ZONE_GROUP);
				if (bgZones != null) {
					bgZones.open();
				}
			}
			catch (Throwable t) {
				bgZones = null;
			}
		}
	}

	/** Prepares the data set. */
	private void prepareExprEntryDataSet() {
		if (entriesExp == null) {
			try {
				entriesExp = genericData.getDataSet(CHP_EXPR_GROUP, CHP_EXPR_GROUP);
				if (entriesExp != null) {
					entriesExp.open();
					wideProbeSetNames = (entriesExp.getHeader().getColumnInfo(0).getColumnType() == DataSetColumnTypes.UnicodeCharColType);
				}
			}
			catch (Throwable t) {
				entriesExp = null;
			}
		}
	}

	/** Prepares the data set for the force call set. */
	private void prepareForceDataSet() {
		if (forceSet == null) {
			try {
				forceSet = genericData.getDataSet(CHP_RESEQ_FORCE_CALL_GROUP, CHP_RESEQ_FORCE_CALL_GROUP);
				if (forceSet != null) {
					forceSet.open();
				}
			}
			catch (Throwable t) {
				forceSet = null;
			}
		}
	}

	/** Prepares the data set. */
	private void prepareGenoEntryDataSet() {
		if (entriesGeno == null) {
			try {
				entriesGeno = genericData.getDataSet(CHP_GENO_GROUP, CHP_GENO_GROUP);
				if (entriesGeno != null) {
					entriesGeno.open();
					wideProbeSetNames = (entriesGeno.getHeader().getColumnInfo(0).getColumnType() == DataSetColumnTypes.UnicodeCharColType);
				}
			}
			catch (Throwable t) {
				entriesGeno = null;
			}
		}
	}

	/** Prepares the data set for the orig call set. */
	private void prepareOrigDataSet() {
		if (origSet == null) {
			try {
				origSet = genericData.getDataSet(CHP_RESEQ_ORIG_CALL_GROUP, CHP_RESEQ_ORIG_CALL_GROUP);
				if (origSet != null) {
					origSet.open();
				}
			}
			catch (Throwable t) {
				origSet = null;
			}
		}
	}

	/** Prepares the data set. */
	private void prepareReseqEntryDataSet() {
		if (entriesReseq == null) {
			try {
				entriesReseq = genericData.getDataSet(CHP_RESEQ_GROUP, CHP_RESEQ_GROUP);
				if (entriesReseq != null) {
					entriesReseq.open();
				}
			}
			catch (Throwable t) {
				entriesReseq = null;
			}
		}
	}

	/** Prepares the data set. */
	private void prepareUnivEntryDataSet() {
		if (entriesUniv == null) {
			try {
				entriesUniv = genericData.getDataSet(CHP_UNIV_GROUP, CHP_UNIV_GROUP);
				if (entriesUniv != null) {
					entriesUniv.open();
				}
			}
			catch (Throwable t) {
				entriesUniv = null;
			}
		}
	}

	public void setBackgroundZoneCnt(int ln) {
		DataSetHeader dpHdr = new DataSetHeader();
		dpHdr.setRowCnt(ln);
		dpHdr.setName(CHP_BG_ZONE_GROUP);

		// center X coord - float
		dpHdr.addFloatColumn(CenterXColName);
		// center Y coord - float
		dpHdr.addFloatColumn(CenterYColName);
		// background - float
		dpHdr.addFloatColumn(BackgroundColName);
		// smoothing factor - float
		dpHdr.addFloatColumn(SmoothFactorColName);

		DataGroupHeader dcHdr = genericData.getHeader().getDataGroup(1);
		dcHdr.addDataSetHdr(dpHdr);
	}

	public void setEntryCount(int ln, int maxln, boolean hasCompData) {
		maxProbeSetName = maxln;
		DataSetHeader dpHdr = new DataSetHeader();
		dpHdr.setRowCnt(ln);
		addColumns(dpHdr, hasCompData);
		DataGroupHeader dcHdr = genericData.getHeader().getDataGroup(0);
		dcHdr.addDataSetHdr(dpHdr);
	}

	/** Sets the name of the CHP file. */
	public void setFilename(String p) {
		genericData.getHeader().setFilename(p);
	}

}
