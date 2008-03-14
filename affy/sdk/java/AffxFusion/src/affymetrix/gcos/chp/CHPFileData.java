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

package affymetrix.gcos.chp;

import java.io.File;
import java.io.FileInputStream;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.channels.FileChannel.MapMode;
import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.utils.IOUtils;
import affymetrix.gcos.FileIO;
import affymetrix.gcos.TagValuePair;
import affymetrix.portability.DataSizes;
import affymetrix.portability.UByte;
import affymetrix.portability.UShort;

/** Provides storage for data in a GCOS CHP file. */
public class CHPFileData {

	/** The CHP file magic number */
	private static final int CHP_FILE_MAGIC_NUMBER = 65;

	/** The max CHP file version the parser can read */
	private static final int CHP_FILE_VERSION_NUMBER = 2;

	/** Identifier to indicate absolute expression analysis results stored. */
	private static final int EXPRESSION_ABSOLUTE_STAT_ANALYSIS = 2;

	/** Identifier to indicate comparison expression analysis results stored. */
	private static final int EXPRESSION_COMPARISON_STAT_ANALYSIS = 3;

	/** Used to convert floating point values stored as ints in older CHP files. */
	private static final int ROUNDFLOAT = 1000;

	/** The file header object */
	private CHPFileHeader header;

	/** The full path of the CHP file */
	private String fileName;

	/** A string to hold an error message associated with a read operation */
	private String strError;

	/** The vector of expression probe set results */
	private List<ExpressionProbeSetResults> exprProbeSetResults = new ArrayList<ExpressionProbeSetResults>();

	/** The vector of universal probe set results */
	private List<UniversalProbeSetResults> univProbeSetResults = new ArrayList<UniversalProbeSetResults>();

	/** The vector of genotype probe set results */
	private List<GenotypeProbeSetResults> genoProbeSetResults = new ArrayList<GenotypeProbeSetResults>();

	/** The resequencing results. */
	private ResequencingResults reseqResults;

	/** Creates a new instance of CHPFileData */
	public CHPFileData() {
		fileName = "";
		clear();
	}

	/**
	 * Opens the file for reading.
	 * 
	 * @param bReadHeaderOnly
	 *          Flag to indicate if the header is to be read only.
	 * @return True if successful.
	 */
	private boolean open(boolean bReadHeaderOnly) {

		// Clear and read the data
		clear();
		if (isXDACompatibleFile()) {
			return readXDAFile(bReadHeaderOnly);
		}
		else {
			return readNonXDAFile(bReadHeaderOnly);
		}
	}

	/**
	 * Read an XDA file.
	 * 
	 * @param bReadHeaderOnly
	 *          Flag to indicate if the header is to be read only.
	 * @return True if successful.
	 */
	private boolean readXDAFile(boolean bReadHeaderOnly) {
		FileInputStream fis = null;
		try {
			// Open the file.
			fis = new FileInputStream(fileName);

			// Read the header.
			if (!readHeaderXDA(fis)) {
				return false;
			}

			// Return if only reading the header.
			if (bReadHeaderOnly == true) {
				return true;
			}

			// Only continue if genotyping or expression
			int assayType = header.getAssayType();
			if ((assayType != CHPFileHeader.EXPRESSION_ASSAY) && (assayType != CHPFileHeader.GENOTYPING_ASSAY)
					&& (assayType != CHPFileHeader.UNIVERSAL_ASSAY) && (assayType != CHPFileHeader.RESEQUENCING_ASSAY)) {
				strError = "The software only supports reading expression, genotyping, tag or resequencing CHP files.";
				return false;
			}

			// Read the probe set data
			if (assayType == CHPFileHeader.EXPRESSION_ASSAY) {
				return readExpressionXDA(fis);
			}
			else if (assayType == CHPFileHeader.GENOTYPING_ASSAY) {
				return readGenotypingXDA(fis);
			}
			else if (assayType == CHPFileHeader.UNIVERSAL_ASSAY) {
				return readUniversalXDA(fis);
			}
			else if (assayType == CHPFileHeader.RESEQUENCING_ASSAY) {
				return readReseqXDA(fis);
			}
			else {
				return false;
			}
		}
		catch (Throwable t) {
			strError = t.getMessage();
			return false;
		}
		finally {
			try {
				if (fis != null) {
					fis.close();
				}
			}
			catch (Exception e) {
			}
		}
	}

	/**
	 * Read the header from an XDA file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return True if successful.
	 */
	private boolean readHeaderXDA(FileInputStream fis) {
		// Read the magic number.
		header = new CHPFileHeader();
		header.setMagic(FileIO.readInt32(fis));

		// Check if new type.
		if (header.getMagic() != CHP_FILE_MAGIC_NUMBER) {
			strError = "The file does not appear to be the correct format.";
			return false;
		}

		// Read the version
		header.setVersion(FileIO.readInt32(fis));

		// Check for version 1 or 2
		if (header.getVersion() > CHP_FILE_VERSION_NUMBER) {
			strError = "Unable to read this version of the CHP file.";
			return false;
		}

		// Get the dimensions of the array
		header.setCols(FileIO.readUInt16(fis));
		header.setRows(FileIO.readUInt16(fis));

		// Number of probe sets.
		header.setNumProbeSets(FileIO.readInt32(fis));
		/* int ival = */FileIO.readInt32(fis); // no qc data extracted.

		// Assay type
		header.setAssayType(FileIO.readInt32(fis));

		// Prog ID.
		header.setProgID(FileIO.readString(fis));

		// Parent cell file.
		header.setParentCellFile(FileIO.readString(fis));

		// Chip type
		header.setChipType(FileIO.readString(fis));

		// Algorithm
		header.setAlgName(FileIO.readString(fis));

		// Algorithm version
		header.setAlgVersion(FileIO.readString(fis));

		// Algorithm parameters.
		int nParams = FileIO.readInt32(fis);
		if (nParams > 0) {
			List<TagValuePair> algParams = new ArrayList<TagValuePair>();
			for (int i = 0; i < nParams; i++) {
				TagValuePair param = new TagValuePair();
				param.setTag(FileIO.readString(fis));
				param.setValue(FileIO.readString(fis));
				algParams.add(param);
			}
			header.setAlgorithmParameters(algParams);
		}

		// Summary parameters
		nParams = FileIO.readInt32(fis);
		if (nParams > 0) {
			List<TagValuePair> chipSummary = new ArrayList<TagValuePair>();
			for (int i = 0; i < nParams; i++) {
				TagValuePair param = new TagValuePair();
				param.setTag(FileIO.readString(fis));
				param.setValue(FileIO.readString(fis));
				chipSummary.add(param);
			}
			header.setSummaryParameters(chipSummary);
		}

		// Background info
		nParams = FileIO.readInt32(fis);
		float sf = FileIO.readFloat(fis);
		if (nParams > 0) {
			BackgroundZoneInfo zones = new BackgroundZoneInfo();
			zones.setSmoothFactor(sf);
			for (int iParam = 0; iParam < nParams; iParam++) {
				BackgroundZoneType zone = new BackgroundZoneType();
				zone.setCenterX(FileIO.readFloat(fis));
				zone.setCenterY(FileIO.readFloat(fis));
				zone.setBackground(FileIO.readFloat(fis));
				zones.addZone(zone);
			}
			header.setBackgroundZoneInfo(zones);
		}
		return true;
	}

	/**
	 * Read the universal data from the XDA file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return True if successful.
	 */
	private boolean readUniversalXDA(FileInputStream fis) {
		// Read each probe set result.
		int n = header.getNumProbeSets();
		univProbeSetResults.clear();
		// int dataSize = FileIO.readInt32(fis);
		for (int i = 0; i < n; i++) {
			UniversalProbeSetResults results = new UniversalProbeSetResults();

			// Read probe set result.
			results.setBackground(FileIO.readFloat(fis));
			univProbeSetResults.add(results);
		}
		return true;
	}

	/**
	 * Read the resequencing data from the XDA file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return True if successful.
	 */
	private boolean readReseqXDA(FileInputStream fis) {
		int dataSize = FileIO.readInt32(fis); // not used.

		// Read the base calls and scores.
		reseqResults = new ResequencingResults();
		dataSize = FileIO.readInt32(fis);
		reseqResults.clearCalledBases();
		reseqResults.clearScores();
		for (int i = 0; i < dataSize; i++) {
			reseqResults.addCalledBase((char)FileIO.readInt8(fis));
		}
		for (int i = 0; i < dataSize; i++) {
			reseqResults.addScore(FileIO.readFloat(fis));
		}

		// Read the force and original calls.
		if (header.getVersion() >= 2) {
			// Read the force calls
			dataSize = FileIO.readInt32(fis);
			reseqResults.clearForceCalls();
			for (int i = 0; i < dataSize; i++) {
				ForceCallType forceCall = new ForceCallType();
				forceCall.setPosition(FileIO.readInt32(fis));
				forceCall.setCall(FileIO.readInt8(fis));
				forceCall.setReason(FileIO.readInt8(fis));
				reseqResults.addForceCall(forceCall);
			}

			// Read the orig calls
			dataSize = FileIO.readInt32(fis);
			reseqResults.clearOrigCalls();
			for (int i = 0; i < dataSize; i++) {
				BaseCallType baseCall = new BaseCallType();
				baseCall.setPosition(FileIO.readInt32(fis));
				baseCall.setCall((char)FileIO.readInt8(fis));
				reseqResults.addOrigCall(baseCall);
			}
		}
		return true;
	}

	/**
	 * Read the expression data from the XDA file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return True if successful.
	 */
	private boolean readExpressionXDA(FileInputStream fis) {
		try {
			// Get the type of analysis
			byte analysisType = FileIO.readInt8(fis); // EXPRESSION_ABSOLUTE_STAT_ANALYSIS
			// or
			// EXPRESSION_COMPARISON_STAT_ANALYSIS
			FileIO.readInt32(fis); // unused
			if ((analysisType != EXPRESSION_ABSOLUTE_STAT_ANALYSIS) && (analysisType != EXPRESSION_COMPARISON_STAT_ANALYSIS)) {
				strError = "The software only supports reading MAS 5 and above expression CHP files.";
				return false;
			}

			// Read each probe set result.
			int n = header.getNumProbeSets();
			exprProbeSetResults.clear();
			for (int i = 0; i < n; i++) {
				ExpressionProbeSetResults results = new ExpressionProbeSetResults();
				// Read the absolute data.
				results.setDetection(new UByte((short)(0xFF & (int)FileIO.readInt8(fis))));
				results.setDetectionPValue(FileIO.readFloat(fis));
				results.setSignal(FileIO.readFloat(fis));
				results.setNumPairs(new UShort(FileIO.readUInt16(fis)));
				results.setNumUsedPairs(new UShort(FileIO.readUInt16(fis)));
				results.setHasCompResults(false);
				// Read the comparison data
				if (analysisType == EXPRESSION_COMPARISON_STAT_ANALYSIS) {
					results.setHasCompResults(true);
					results.setChange(new UByte((short)(0xFF & (int)FileIO.readInt8(fis))));
					results.setChangePValue(FileIO.readFloat(fis));
					results.setSignalLogRatio(FileIO.readFloat(fis));
					results.setSignalLogRatioLow(FileIO.readFloat(fis));
					results.setSignalLogRatioHigh(FileIO.readFloat(fis));
					results.setNumCommonPairs(new UShort(FileIO.readUInt16(fis)));
				}
				exprProbeSetResults.add(results);
			}
		}
		catch (UnsignedOutOfLimitsException e) {
			e.printStackTrace();
		}
		return true;
	}

	/**
	 * Read the genotyping data from the XDA file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return True if successful.
	 */
	private boolean readGenotypingXDA(FileInputStream fis) {
		final int DM_ALG_RESULT_SIZE = 21;

		// Read each probe set result.
		int n = header.getNumProbeSets();
		genoProbeSetResults.clear();
		int dataSize = FileIO.readInt32(fis);
		for (int iset = 0; iset < n; iset++) {
			GenotypeProbeSetResults results = new GenotypeProbeSetResults();

			// Read probe set result.
			results.setAlleleCall(FileIO.readInt8(fis));

			results.setConfidence(FileIO.readFloat(fis));

			results.setRAS1(FileIO.readFloat(fis));
			results.setPValue_AA(results.getRAS1());

			results.setRAS2(FileIO.readFloat(fis));
			results.setPValue_AB(results.getRAS2());

			if (dataSize == DM_ALG_RESULT_SIZE) {
				results.setPValue_BB(FileIO.readFloat(fis));
				results.setPValue_NoCall(FileIO.readFloat(fis));
			}
			genoProbeSetResults.set(iset, results);
		}
		return true;
	}

	/**
	 * Read the non-XDA file header.
	 * 
	 * @param fis
	 *          The file input stream.
	 * @return True if successful.
	 */
	private boolean readHeaderNonXDA(FileInputStream fis, List<Integer> tagCells) {
		// Read the string that defines the CHP file (older format).
		String appName = "GeneChip Sequence File";
		String strMagic = FileIO.readFixedString(fis, appName.length());
		if (strMagic.compareTo(appName) != 0) {
			strError = "The file does not appear to be the correct format.";
			return false;
		}

		// Read the version number
		header = new CHPFileHeader();
		header.setVersion(FileIO.readInt32(fis));

		// Read algorithm type string.
		header.setAlgName(FileIO.readString(fis));
		header.setAlgVersion(FileIO.readString(fis));

		// Read parameters.
		String strParam = FileIO.readString(fis);
		header.setAlgorithmParameters(parseString(strParam, " ", "="));

		// Read summary
		strParam = FileIO.readString(fis);
		header.setSummaryParameters(parseString(strParam, " ", "="));

		// rows
		header.setRows(FileIO.readInt32(fis));

		// cols
		header.setCols(FileIO.readInt32(fis));

		// #probe sets
		header.setNumProbeSets(FileIO.readInt32(fis));

		// the maximum probe set number in the array
		int maxvalue = FileIO.readInt32(fis);

		// #qc probe sets
		/* int nqc = */FileIO.readInt32(fis);

		// probe set numbers
		int ival;
		int nsets = header.getNumProbeSets();
		for (int j = 0; j < nsets; j++) {
			ival = FileIO.readInt32(fis);
		}

		// #probe pairs in each probe set.
		for (int j = 0; j < maxvalue; j++) {
			ival = FileIO.readInt32(fis);
		}

		// type of probe set
		for (int j = 0; j < maxvalue; j++) {
			ival = FileIO.readInt32(fis);
			if (j == 0) {
				switch (ival) {
				case 1:
					header.setAssayType(CHPFileHeader.RESEQUENCING_ASSAY);
					break;

				case 2:
					header.setAssayType(CHPFileHeader.GENOTYPING_ASSAY);
					break;

				case 3:
					header.setAssayType(CHPFileHeader.EXPRESSION_ASSAY);
					break;

				case 7:
					header.setAssayType(CHPFileHeader.UNIVERSAL_ASSAY);
					break;

				default:
					header.setAssayType(CHPFileHeader.UNKNOWN_ASSAY);
					break;
				}
			}
		}

		// This must be a resequencing design if there are no probe sets.
		if (nsets == 0) {
			header.setAssayType(CHPFileHeader.RESEQUENCING_ASSAY);
		}

		// Check for valid versions.
		if (((header.getVersion() < 12) && ((header.getAssayType() == CHPFileHeader.EXPRESSION_ASSAY) || (header
				.getAssayType() == CHPFileHeader.GENOTYPING_ASSAY)))
				|| ((header.getVersion() < 10) && (header.getAssayType() == CHPFileHeader.UNIVERSAL_ASSAY))) {
			strError = "This version of the CHP file is not supported by the parser.";
			return false;
		}

		// # probe probes per element (probe pair, probe quartet, etc).
		boolean readTagCells = false;
		if (header.getAssayType() == CHPFileHeader.UNIVERSAL_ASSAY) {
			readTagCells = true;
		}
		for (int j = 0; j < nsets; j++) {
			ival = FileIO.readInt32(fis);
			if (readTagCells == true) {
				tagCells.add(j, ival);
			}
		}

		// Get the chip type
		header.setChipType(FileIO.readFixedString(fis, 256));

		// The parent CEL file.
		header.setParentCellFile(FileIO.readFixedString(fis, 256));

		// The prog ID.
		header.setProgID(FileIO.readString(fis));

		return true;
	}

	/**
	 * Reads the expression data from a non-XDA file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return True if successful.
	 */
	private boolean readExpressionNonXDA(FileInputStream fis) {
		try {
			// Read each probe set result.
			int ival;
			// float fval;
			// byte bval;
			int nsets = header.getNumProbeSets();
			exprProbeSetResults.clear();
			for (int iset = 0; iset < nsets; iset++) {
				ExpressionProbeSetResults results = new ExpressionProbeSetResults();

				results.setNumPairs(new UShort(FileIO.readInt32(fis)));
				results.setNumUsedPairs(new UShort(FileIO.readInt32(fis)));
				if (header.getVersion() <= 12) {
					ival = FileIO.readInt32(fis);
				}

				ival = FileIO.readInt32(fis); // unused
				if (header.getVersion() == 12) {
					ival = FileIO.readInt32(fis); // unused
					ival = FileIO.readInt32(fis); // unused
					ival = FileIO.readInt32(fis); // unused
				}

				results.setDetectionPValue(FileIO.readFloat(fis));

				if (header.getVersion() == 12) {
					/* fval = */FileIO.readFloat(fis); // unused
				}

				results.setSignal(FileIO.readFloat(fis));
				results.setDetection(new UByte((short)FileIO.readInt32(fis)));

				// unused
				int npairs = results.getNumPairs().toInt();
				for (int ip = 0; ip < npairs; ++ip) {
					/* fval = */FileIO.readFloat(fis);
					ival = FileIO.readInt32(fis);
					if (header.getVersion() == 12) {
						ival = FileIO.readInt32(fis);
						ival = FileIO.readInt32(fis);
						/* fval = */FileIO.readFloat(fis);
						/* fval = */FileIO.readFloat(fis);
						ival = FileIO.readInt32(fis);
						/* bval = */FileIO.readInt8(fis);
						/* bval = */FileIO.readInt8(fis);
					}
					else {
						ival = FileIO.readUInt16(fis);
						ival = FileIO.readUInt16(fis);
					}
					if (header.getVersion() == 12) {
						ival = FileIO.readInt32(fis);
						ival = FileIO.readInt32(fis);
						/* fval = */FileIO.readFloat(fis);
						/* fval = */FileIO.readFloat(fis);
						ival = FileIO.readInt32(fis);
						/* bval = */FileIO.readInt8(fis);
						/* bval = */FileIO.readInt8(fis);
					}
					else {
						ival = FileIO.readUInt16(fis);
						ival = FileIO.readUInt16(fis);
					}
				}

				ival = FileIO.readInt32(fis);
				results.setHasCompResults(ival == 1 ? true : false);

				if (results.getHasCompResults()) {
					results.setNumCommonPairs(new UShort(FileIO.readInt32(fis)));

					if (header.getVersion() == 12) {
						ival = FileIO.readInt32(fis); // unused
						ival = FileIO.readInt32(fis); // unused
						ival = FileIO.readInt32(fis); // unused
					}

					results.setChange(new UByte((short)FileIO.readInt32(fis)));

					/* bval = */FileIO.readInt8(fis); // unused
					if (header.getVersion() == 12) {
						/* bval = */FileIO.readInt8(fis); // unused
						ival = FileIO.readInt32(fis); // unused
						ival = FileIO.readInt32(fis); // unused
					}

					results.setSignalLogRatioHigh((float)FileIO.readInt32(fis) / ROUNDFLOAT);

					ival = FileIO.readInt32(fis); // unused

					if (header.getVersion() == 12) {
						ival = FileIO.readInt32(fis); // unused
					}

					results.setSignalLogRatio((float)FileIO.readInt32(fis) / ROUNDFLOAT);

					if (header.getVersion() == 12) {
						ival = FileIO.readInt32(fis); // unused
					}

					results.setSignalLogRatioLow((float)FileIO.readInt32(fis) / ROUNDFLOAT);

					if (header.getVersion() == 12) {
						results.setChangePValue((float)FileIO.readInt32(fis) / ROUNDFLOAT);
					}
					else {
						results.setChangePValue(FileIO.readFloat(fis));
					}
				}
				exprProbeSetResults.set(iset, results);
			}
		}
		catch (UnsignedOutOfLimitsException e) {
			e.printStackTrace();
		}
		return true;
	}

	/**
	 * Reads the genotyping data from a non-XDA file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return True if successful.
	 */
	private boolean readGenotypingNonXDA(FileInputStream fis) {
		int nsets = header.getNumProbeSets();
		genoProbeSetResults.clear();
		// int ival;
		// float fval;
		byte bval;
		// String sval;
		for (int iset = 0; iset < nsets; iset++) {
			GenotypeProbeSetResults results = new GenotypeProbeSetResults();

			// Unused data
			int ngroups = FileIO.readInt32(fis);
			for (int ig = 0; ig < ngroups; ig++) {
				/* ival = */FileIO.readInt32(fis);
				/* sval = */FileIO.readString(fis);
				bval = FileIO.readInt8(fis);
				/* ival = */FileIO.readInt32(fis);
				/* ival = */FileIO.readInt32(fis);
				/* ival = */FileIO.readInt32(fis);
			}
			bval = FileIO.readInt8(fis);
			if (bval == 1) {
				/* ival = */FileIO.readInt32(fis);
				/* sval = */FileIO.readString(fis);
				/* sval = */FileIO.readString(fis);
				/* sval = */FileIO.readString(fis);
				/* ival = */FileIO.readInt32(fis);
				/* ival = */FileIO.readInt32(fis);

				// The call
				results.setAlleleCall(FileIO.readInt8(fis));

				// The confidence
				if (header.getVersion() == 12) {
					results.setConfidence((float)FileIO.readInt32(fis) / ROUNDFLOAT);
				}
				else {
					results.setConfidence(FileIO.readFloat(fis));
				}

				// unused
				/* fval = */FileIO.readFloat(fis);
				/* fval = */FileIO.readFloat(fis);
				/* fval = */FileIO.readFloat(fis);

				// RAS 1 and 2
				results.setRAS1(FileIO.readFloat(fis));
				results.setRAS2(FileIO.readFloat(fis));
			}
			else {
				results.setConfidence(0.0f);
				results.setRAS1(0.0f);
				results.setRAS2(0.0f);
				results.setAlleleCall(GenotypeProbeSetResults.ALLELE_NO_CALL);
			}

			// 100K results are not stored in this version
			results.setPValue_AA(0.0f);
			results.setPValue_AB(0.0f);
			results.setPValue_BB(0.0f);
			results.setPValue_NoCall(0.0f);

			// unused
			/* sval = */FileIO.readString(fis);
			/* sval = */FileIO.readString(fis);
			int np = FileIO.readInt32(fis);
			for (int ip = 0; ip < np; ++ip) {
				/* ival = */FileIO.readInt32(fis);
				if (header.getVersion() == 12) {
					/* ival = */FileIO.readInt32(fis);
					/* ival = */FileIO.readInt32(fis);
					/* ival = */FileIO.readInt32(fis);
					/* ival = */FileIO.readInt32(fis);
					/* ival = */FileIO.readInt32(fis);
					bval = FileIO.readInt8(fis);
					bval = FileIO.readInt8(fis);
				}
				else {
					/* ival = */FileIO.readUInt16(fis);
					/* ival = */FileIO.readUInt16(fis);
				}
				if (header.getVersion() == 12) {
					/* ival = */FileIO.readInt32(fis);
					/* ival = */FileIO.readInt32(fis);
					/* ival = */FileIO.readInt32(fis);
					/* ival = */FileIO.readInt32(fis);
					/* ival = */FileIO.readInt32(fis);
					bval = FileIO.readInt8(fis);
					bval = FileIO.readInt8(fis);
				}
				else {
					/* ival = */FileIO.readUInt16(fis);
					/* ival = */FileIO.readUInt16(fis);
				}
			}
			genoProbeSetResults.set(iset, results);
		}
		return true;
	}

	/**
	 * Reads the universal data from a non-XDA file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @param tagCells
	 *          The number of cells per atom.
	 * @return True if successful.
	 */
	private boolean readUniversalNonXDA(FileInputStream fis, List<Integer> tagCells) {
		int nsets = header.getNumProbeSets();
		univProbeSetResults.clear();
		// byte bval;
		int ival;

		// Read each probe set result.
		for (int i = 0; i < nsets; i++) {
			UniversalProbeSetResults results = new UniversalProbeSetResults();

			// unused (wildtype) length(int), string(len)
			ival = FileIO.readInt32(fis);
			for (int n = 0; n < ival; n++) {
				FileIO.readInt8(fis);
			}
			// unused (called) length(int), string(len)
			ival = FileIO.readInt32(fis);
			for (int n = 0; n < ival; n++) {
				FileIO.readInt8(fis);
			}

			// unused (#atoms(int))
			int nAtoms = FileIO.readInt32(fis);
			for (int iatom = 0; iatom < nAtoms; iatom++) {
				ival = FileIO.readInt32(fis);
				if (iatom == 0) {
					results.setBackground((float)ival / ROUNDFLOAT);
				}

				// Ignore the probe level data.
				int ncells = tagCells.get(i);
				for (int n = 0; n < ncells; n++) {
					if (header.getVersion() <= 12) {
						ival = FileIO.readInt32(fis); // unused X coordinate
						ival = FileIO.readInt32(fis); // unused Y coordinate
					}
					else {
						ival = FileIO.readUInt16(fis); // unused X
						// coordinate
						ival = FileIO.readUInt16(fis); // unused Y
						// coordinate
					}
					if (header.getVersion() <= 12) {
						ival = FileIO.readInt32(fis); // unused intensity
						ival = FileIO.readInt32(fis); // unused stdv
						ival = FileIO.readInt32(fis); // unused pixel count
						if (header.getVersion() >= 8) {
							FileIO.readInt8(fis); // unused mask
							FileIO.readInt8(fis); // unused
							// outlier
						}
					}
				}
			}
			univProbeSetResults.add(results);
		}
		return true;
	}

	/**
	 * Reads the resequencing data from a non-XDA file.
	 * 
	 * @param fis
	 *          The input file stream.
	 * @return True if successful.
	 */
	private boolean readReseqNonXDA(FileInputStream fis) {
		int ival;
		// Read the base calls.
		reseqResults = new ResequencingResults();
		int dataSize = FileIO.readInt32(fis);
		reseqResults.clearCalledBases();
		reseqResults.clearScores();
		for (int i = 0; i < dataSize; i++) {
			reseqResults.addCalledBase((char)FileIO.readInt8(fis));
		}
		// Ignore the edit section - there were none with the non-XDA file.
		for (int i = 0; i < dataSize; i++) {
			FileIO.readInt8(fis);
		}
		// Ignore the confidence section - there were none with the non-XDA file.
		for (int i = 0; i < dataSize; i++) {
			FileIO.readInt8(fis);
		}
		// Ignore the unit index section - there were none with the non-XDA file.
		for (int i = 0; i < dataSize; i++) {
			ival = FileIO.readInt16(fis);
		}
		ival = FileIO.readInt32(fis);
		if (ival > 0) {
			FileIO.readFixedString(fis, ival);
			ival = FileIO.readInt32(fis);
			if (ival > 0) {
				FileIO.readFixedString(fis, ival);
			}
		}
		for (int i = 0; i < dataSize; i++) {
			reseqResults.addScore(FileIO.readFloat(fis));
		}
		return true;
	}

	/**
	 * Read an non-XDA file.
	 * 
	 * @param bReadHeaderOnly
	 *          Flag to indicate if the header is to be read only.
	 * @return True if successful.
	 */
	private boolean readNonXDAFile(boolean bReadHeaderOnly) {
		FileInputStream fis = null;
		try {
			// Open the file.
			fis = new FileInputStream(fileName);

			// Read the header.
			List<Integer> tagCells = new ArrayList<Integer>();
			if (!readHeaderNonXDA(fis, tagCells)) {
				return false;
			}
			// Return if only reading the header.
			if (bReadHeaderOnly == true) {
				return true;
			}
			// Only continue if genotyping or expression
			int assayType = header.getAssayType();
			if ((assayType != CHPFileHeader.EXPRESSION_ASSAY) && (assayType != CHPFileHeader.GENOTYPING_ASSAY)
					&& (assayType != CHPFileHeader.UNIVERSAL_ASSAY) && (assayType != CHPFileHeader.RESEQUENCING_ASSAY)) {
				strError = "The software only supports reading expression, genotyping, tag or resequencing CHP files.";
				return false;
			}

			// Read the probe set data
			if (assayType == CHPFileHeader.EXPRESSION_ASSAY) {
				return readExpressionNonXDA(fis);

			}
			else if (assayType == CHPFileHeader.GENOTYPING_ASSAY) {
				return readGenotypingNonXDA(fis);
			}
			else if (assayType == CHPFileHeader.UNIVERSAL_ASSAY) {
				return readUniversalNonXDA(fis, tagCells);
			}
			else if (assayType == CHPFileHeader.RESEQUENCING_ASSAY) {
				return readReseqNonXDA(fis);
			}
			else {
				return false;
			}
		}
		catch (Throwable t) {
			strError = t.getMessage();
			return false;
		}
		finally {
			try {
				if (fis != null) {
					fis.close();
				}
			}
			catch (Exception e) {
			}
		}
	}

	/**
	 * Accessors to header.
	 * 
	 * @return The header data object
	 */
	public CHPFileHeader getHeader() {
		return header;
	}

	/**
	 * Returns the expression probe set result
	 * 
	 * @param index
	 *          The index to the result object of interest.
	 * @return The expression result.
	 */
	public ExpressionProbeSetResults getExpressionResults(int index) {
		if ((index < header.getNumProbeSets()) && (header.getAssayType() == CHPFileHeader.EXPRESSION_ASSAY)) {
			return exprProbeSetResults.get(index);
		}
		return null;
	}

	/**
	 * Returns the genotyping probe set result
	 * 
	 * @param index
	 *          The index to the result object of interest.
	 * @return The genotyping result.
	 */
	public GenotypeProbeSetResults getGenotypingResults(int index) {
		if ((index < header.getNumProbeSets()) && (header.getAssayType() == CHPFileHeader.GENOTYPING_ASSAY)) {
			return genoProbeSetResults.get(index);
		}
		return null;
	}

	/**
	 * Returns the universal (tag array) probe set result
	 * 
	 * @param index
	 *          The index to the result object of interest.
	 * @return The universal result.
	 */
	public UniversalProbeSetResults getUniversalResults(int index) {
		if ((index < header.getNumProbeSets()) && (header.getAssayType() == CHPFileHeader.UNIVERSAL_ASSAY)) {
			return univProbeSetResults.get(index);
		}
		return null;
	}

	/**
	 * Returns the resequencing results.
	 * 
	 * @return The resequencing results.
	 */
	public ResequencingResults getResequencingResults() {
		if (header.getAssayType() == CHPFileHeader.RESEQUENCING_ASSAY) {
			return reseqResults;
		}
		return null;
	}

	/**
	 * Error string when the read functions fail.
	 * 
	 * @return A string message describing a read error
	 */
	public String getError() {
		return strError;
	}

	/**
	 * Reads the entire file.
	 * 
	 * @return True if successful.
	 */
	public boolean read() {
		// Read the header, close if failed.
		if (!open(false)) {
			clear();
			return false;
		}
		return true;
	}

	/**
	 * Reads the header of the CHP file
	 * 
	 * @return True if successful
	 */
	public boolean readHeader() {
		// Read the header, close if failed.
		if (open(true) == false) {
			clear();
			return false;
		}
		return true;
	}

	/**
	 * Determines if the file specified by the FileName property exists.
	 * 
	 * @return True if the file exists.
	 */
	public boolean exists() {
		return new File(fileName).exists();
	}

	/**
	 * Determines if the CHP file specified by the FileName property is an XDA format file
	 * 
	 * @return True if the file is an XDA file
	 */
	public boolean isXDACompatibleFile() {

		FileInputStream fis = null;
		// Open the file.
		try {
			fis = new FileInputStream(fileName);
			FileChannel fc = fis.getChannel();
			MappedByteBuffer mbb = fc.map(MapMode.READ_ONLY, 0, DataSizes.INT_SIZE);
			mbb.order(ByteOrder.LITTLE_ENDIAN);
			int magic = mbb.getInt();
			mbb = null;
			return (magic == CHPFileData.CHP_FILE_MAGIC_NUMBER);
		}
		catch (Throwable t) {
			return false;
		}
		finally {
			try {
				if (fis != null) {
					fis.close();
				}
			}
			catch (Exception e) {
			}
		}
	}

	/**
	 * Sets the file name.
	 * 
	 * @param name
	 *          The full path to the CHP file
	 */
	public void setFileName(String name) {
		fileName = name;
	}

	/**
	 * Gets the file name.
	 * 
	 * @return The full path to the CHP file.
	 */
	public String getFileName() {
		return fileName;
	}

	/** Deallocates any memory used by the class object */
	public void clear() {
		header = null;
		strError = IOUtils.EMPTY;
		exprProbeSetResults.clear();
		genoProbeSetResults.clear();
		univProbeSetResults.clear();
		reseqResults = null;
	}

	/**
	 * Parses a string into tag/value parameters.
	 * 
	 * @param sSource
	 *          The parameter string
	 * @sDelimiter1 The delimiter between the tag and value
	 * @sDelimiter2 The delimiter between the parameters.
	 */
	private List<TagValuePair> parseString(String source, String delim1, String delim2) {

		if (source == null) {
			return null;
		}

		// Parameters are like one of the following (no value, spaces in values,
		// simple tag/values)
		// BF= NF=1
		// BF=c:\p f\file.cel NF=1

		List<TagValuePair> tagValueParams = new ArrayList<TagValuePair>();
		String[] params = source.split(delim1);
		int nParams = 0;
		int index;
		for (int i = 0; i < params.length; i++) {
			index = params[i].indexOf(delim2);

			// No value, just the tag.
			if (index == params[i].length() - 1) {
				continue;
			}
			// Tag and value
			else if (index >= 0) {
				String[] tagvalue = params[i].split(delim2);
				TagValuePair param = new TagValuePair();
				param.setTag(tagvalue[0]);
				param.setValue(tagvalue[1]);
				tagValueParams.add(param);
				++nParams;
			}

			// No tag, just value.
			else if (index == -1) {
				TagValuePair param = tagValueParams.get(nParams - 1);
				param.setValue(param.getValue() + " " + params[i]);
			}
		}
		return tagValueParams;
	}

}
