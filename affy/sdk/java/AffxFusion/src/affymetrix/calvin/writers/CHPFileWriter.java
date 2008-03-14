/////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
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
package affymetrix.calvin.writers;

import java.io.IOException;
import java.util.Iterator;

import affymetrix.calvin.data.CHPBackgroundZone;
import affymetrix.calvin.data.CHPData;
import affymetrix.calvin.data.CHPExpressionEntry;
import affymetrix.calvin.data.CHPGenotypeEntry;
import affymetrix.calvin.data.CHPReseqEntry;
import affymetrix.calvin.data.CHPReseqForceCall;
import affymetrix.calvin.data.CHPReseqOrigCall;
import affymetrix.calvin.data.CHPUniversalEntry;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.portability.UInt;

/**
 * The Class CHPFileWriter.
 */
public class CHPFileWriter extends CHPFileWriterBase {

	/** The maximum length of a probe set name. */
	int maxProbeSetName;

	/** The data set writer. */
	private DataSetWriter dataSetWriter;

	/** The position of the data entry. */
	private UInt entryPos = UInt.ZERO;

	/** The position of the bg zone entry. */
	private UInt bgZonePos = UInt.ZERO;

	/** The position of the force entry. */
	private UInt forcePos = UInt.ZERO;

	/** The position of the orig entry. */
	private UInt origPos = UInt.ZERO;

	/**
	 * Instantiates a new CHP file writer.
	 * 
	 * @param p
	 *          the CHP data
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 * @throws UnsignedOutOfLimitsException
	 *           the unsigned out of limits exception
	 */
	CHPFileWriter(CHPData p) throws IOException, UnsignedOutOfLimitsException {
		super(p.getFileHeader());
		maxProbeSetName = p.getMaxProbeSetName();
		writeHeaders();
	}

	/**
	 * Write file headers.
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 * @throws UnsignedOutOfLimitsException
	 *           the unsigned out of limits exception
	 */
	private void writeHeaders() throws IOException, UnsignedOutOfLimitsException {
		writer.writeHeader();
		Iterator<DataGroupWriter> dataGrpIt = writer.getDataGroupWriterIterator();
		while (dataGrpIt.hasNext()) {
			DataGroupWriter dataGroupWriter = dataGrpIt.next();
			dataGroupWriter.writeHeader();
			Iterator<DataSetWriter> dataSetIt = dataGroupWriter.getDataSetWriterIterator();
			while (dataSetIt.hasNext()) {
				dataSetWriter = dataSetIt.next();
				dataSetWriter.writeHeader();
				setFilePositions();
			}
			dataGroupWriter.updateNextDataGroupPos();
		}
	}

	/**
	 * Seek to data set.
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	public void seekToDataSet() throws IOException {
		writer.seekFromBeginPos(entryPos);
	}

	/**
	 * Write expression entry.
	 * 
	 * @param p
	 *          The entry to write.
	 * 
	 * @throws Exception
	 *           the exception
	 */
	public void writeExpressionEntry(CHPExpressionEntry p) throws /* IOException */Exception {
		dataSetWriter.write8Bit(p.getProbeSetName(), maxProbeSetName);
		dataSetWriter.write(p.getDetection());
		dataSetWriter.write(p.getDetectionPValue());
		dataSetWriter.write(p.getSignal());
		dataSetWriter.write(p.getNumPairs());
		dataSetWriter.write(p.getNumPairsUsed());
		if (p.getHasComparisonData()) {
			dataSetWriter.write(p.getChange());
			dataSetWriter.write(p.getChangePValue());
			dataSetWriter.write(p.getSigLogRatio());
			dataSetWriter.write(p.getSigLogRatioLo());
			dataSetWriter.write(p.getSigLogRatioHi());
			dataSetWriter.write(p.getCommonPairs());
		}
	}

	/**
	 * Write genotype entry.
	 * 
	 * @param p
	 *          The entry to write.
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	public void writeGenotypeEntry(CHPGenotypeEntry p) throws IOException {
		dataSetWriter.write8Bit(p.getProbeSetName(), maxProbeSetName);
		dataSetWriter.write(p.getCall());
		dataSetWriter.write(p.getConfidence());
		dataSetWriter.write(p.getRAS1());
		dataSetWriter.write(p.getRAS2());
		dataSetWriter.write(p.getAACall());
		dataSetWriter.write(p.getABCall());
		dataSetWriter.write(p.getBBCall());
		dataSetWriter.write(p.getNoCall());
	}

	/**
	 * Write universal entry.
	 * 
	 * @param p
	 *          The entry to write.
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	public void writeUniversalEntry(CHPUniversalEntry p) throws IOException {
		dataSetWriter.write(p.getBackground());
	}

	/**
	 * Write reseq entry.
	 * 
	 * @param p
	 *          The entry to write.
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	public void writeReseqEntry(CHPReseqEntry p) throws IOException {
		dataSetWriter.write(p.getCall());
		dataSetWriter.write(p.getScore());
	}

	/**
	 * Seek to bg set.
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	public void seekToBgSet() throws IOException {
		writer.seekFromBeginPos(bgZonePos);
	}

	/**
	 * Write background zone.
	 * 
	 * @param zone
	 *          The zone data to write.
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	public void writeBackgroundZone(CHPBackgroundZone zone) throws IOException {
		dataSetWriter.write(zone.getCenterX());
		dataSetWriter.write(zone.getCenterY());
		dataSetWriter.write(zone.getBackground());
		dataSetWriter.write(zone.getSmoothFactor());
	}

	/**
	 * Seek to force set.
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	public void seekToForceSet() throws IOException {
		writer.seekFromBeginPos(forcePos);
	}

	/**
	 * Write force call.
	 * 
	 * @param force
	 *          the force call
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	public void writeForceCall(CHPReseqForceCall force) throws IOException {
		dataSetWriter.write(force.getPosition());
		dataSetWriter.write(force.getCall());
		dataSetWriter.write(force.getReason());
	}

	/**
	 * Seek to orig call data set.
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	public void seekToOrigCallSet() throws IOException {
		writer.seekFromBeginPos(origPos);
	}

	/**
	 * Write orig call.
	 * 
	 * @param orig
	 *          the orig call
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	public void writeOrigCall(CHPReseqOrigCall orig) throws IOException {
		dataSetWriter.write(orig.getPosition());
		dataSetWriter.write(orig.getCall());
	}

	/**
	 * Updates the file positions after a data set is added to the file.
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 * @throws UnsignedOutOfLimitsException
	 *           the unsigned out of limits exception
	 */
	private void setFilePositions() throws IOException, UnsignedOutOfLimitsException {
		String name = dataSetWriter.getName();
		int dataSetSz = dataSetWriter.getSize();
		if (name == CHPData.CHP_BG_ZONE_GROUP) {
			bgZonePos = writer.getFilePos();
		}
		else if (name == CHPData.CHP_RESEQ_FORCE_CALL_GROUP) {
			forcePos = writer.getFilePos();
		}
		else if (name == CHPData.CHP_RESEQ_ORIG_CALL_GROUP) {
			origPos = writer.getFilePos();
		}
		else {
			entryPos = writer.getFilePos();
		}

		writer.seekFromCurrentPos(new UInt(0xFFFFFFFFL & (long)(dataSetSz + 1)));
		dataSetWriter.updateNextDataSetOffset();
	}
}
