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
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import affymetrix.calvin.data.CHPMultiDataData;
import affymetrix.calvin.data.DataSetInfo;
import affymetrix.calvin.data.ProbeSetMultiDataCopyNumberData;
import affymetrix.calvin.data.ProbeSetMultiDataCytoRegionData;
import affymetrix.calvin.data.ProbeSetMultiDataExpressionData;
import affymetrix.calvin.data.ProbeSetMultiDataGenotypeData;
import affymetrix.calvin.data.CHPMultiDataData.MultiDataType;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.portability.UInt;

/**
 * The Class CHPMultiDataFileWriter.
 */
public class CHPMultiDataFileWriter extends CHPFileWriterBase {

	/** The data set writer. */
	private DataSetWriter dataSetWriter;

	/** The file position of the entry for each data set. */
	private Map<MultiDataType, UInt> entryPos = new HashMap<MultiDataType, UInt>();

	/** A map of data type to data set index. */
	private Map<MultiDataType, Integer> dataTypeToIndex = new HashMap<MultiDataType, Integer>();

	/** A map of data set index to data type. */
	private Map<Integer, MultiDataType> indexToDataType = new HashMap<Integer, MultiDataType>();

	/** The maximum probe set name. */
	private int maxProbeSetName;

	/** The data type being written. */
	private MultiDataType currentDataType;

	/** The data. */
	private CHPMultiDataData data;

	/**
	 * Instantiates a new CHP multi data file writer.
	 * 
	 * @param p
	 *          multi data
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 * @throws UnsignedOutOfLimitsException
	 *           the unsigned out of limits exception
	 */
	public CHPMultiDataFileWriter(CHPMultiDataData p) throws IOException, UnsignedOutOfLimitsException {
		super(p.getFileHeader());
		data = p;
		Map<MultiDataType, DataSetInfo> info = p.getDataSetInfo();
		Iterator<Entry<MultiDataType, DataSetInfo>> it = info.entrySet().iterator();
		while (it.hasNext()) {
			Map.Entry<MultiDataType, DataSetInfo> entry = it.next();
			dataTypeToIndex.put(entry.getKey(), entry.getValue().getDataSetIndex());
			indexToDataType.put(entry.getValue().getDataSetIndex(), entry.getKey());
		}
		writeHeaders();
	}

	/**
	 * Seek to data set.
	 * 
	 * @param dataType
	 *          the data type
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	public void seekToDataSet(MultiDataType dataType) throws IOException {
		DataGroupWriter dataGroupWriter = writer.getDataGroupWriter(0);
		writer.seekFromBeginPos(entryPos.get(dataType));
		dataSetWriter = dataGroupWriter.getDataSetWriter(dataTypeToIndex.get(dataType));
		maxProbeSetName = data.getMaxProbeSetName(dataType);
		currentDataType = dataType;
	}

	/**
	 * Write entry.
	 * 
	 * @param p
	 *          the ProbeSetMultiDataGenotypeData
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 * @throws UnsignedOutOfLimitsException
	 *           the unsigned out of limits exception
	 */
	public void writeEntry(ProbeSetMultiDataGenotypeData p) throws IOException, UnsignedOutOfLimitsException {

		writer.seekFromBeginPos(entryPos.get(currentDataType));
		dataSetWriter.write8Bit(p.getName(), maxProbeSetName);
		dataSetWriter.write(p.getCall());
		dataSetWriter.write(p.getConfidence());
		writeMetrics(p.getMetrics());
		entryPos.put(currentDataType, writer.getFilePos());
	}

	/**
	 * Write entry.
	 * 
	 * @param p
	 *          the ProbeSetMultiDataCopyNumberData
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 * @throws UnsignedOutOfLimitsException
	 *           the unsigned out of limits exception
	 */
	public void writeEntry(ProbeSetMultiDataCopyNumberData p) throws IOException, UnsignedOutOfLimitsException {
		writer.seekFromBeginPos(entryPos.get(currentDataType));
		dataSetWriter.write8Bit(p.getName(), maxProbeSetName);
		dataSetWriter.write(p.getChr());
		dataSetWriter.write(p.getPosition());
		writeMetrics(p.getMetrics());
		entryPos.put(currentDataType, writer.getFilePos());
	}

	public void writeEntry(ProbeSetMultiDataCytoRegionData p) throws IOException, UnsignedOutOfLimitsException {
		writer.seekFromBeginPos(entryPos.get(currentDataType));
		dataSetWriter.write8Bit(p.getName(), maxProbeSetName);
		dataSetWriter.write(p.getChr());
		dataSetWriter.write(p.getStartPosition());
		dataSetWriter.write(p.getStopPosition());
		dataSetWriter.write(p.getCall());
		dataSetWriter.write(p.getConfidence());
		writeMetrics(p.getMetrics());
		entryPos.put(currentDataType, writer.getFilePos());
	}

	/**
	 * Write entry.
	 * 
	 * @param p
	 *          the multi data expression
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 * @throws UnsignedOutOfLimitsException
	 *           the unsigned out of limits exception
	 */
	public void writeEntry(ProbeSetMultiDataExpressionData p) throws IOException, UnsignedOutOfLimitsException {
		writer.seekFromBeginPos(entryPos.get(currentDataType));
		dataSetWriter.write8Bit(p.getName(), maxProbeSetName);
		dataSetWriter.write(p.getQuantification());
		writeMetrics(p.getMetrics());
		entryPos.put(currentDataType, writer.getFilePos());
	}

	/**
	 * Write metrics.
	 * 
	 * @param metrics
	 *          the metrics
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	private void writeMetrics(List<ParameterNameValue> metrics) throws IOException, UnsignedOutOfLimitsException {
		int ncols = metrics.size();
		for (int icol = 0; icol < ncols; icol++) {
			ParameterNameValue nv = metrics.get(icol);
			switch (nv.getParameterType()) {
			case Int8Type:
				dataSetWriter.write(nv.getValueInt8());
				break;

			case UInt8Type:
				dataSetWriter.write(nv.getValueUInt8());
				break;

			case Int16Type:
				dataSetWriter.write(nv.getValueInt16());
				break;

			case UInt16Type:
				dataSetWriter.write(nv.getValueUInt16());
				break;

			case Int32Type:
				dataSetWriter.write(nv.getValueInt32());
				break;

			case UInt32Type:
				dataSetWriter.write(nv.getValueUInt32());
				break;

			case FloatType:
				dataSetWriter.write(nv.getValueFloat());
				break;

			case AsciiType:
				dataSetWriter.write8Bit(nv.getValueAscii(), data.getMetricColumnLength(currentDataType, icol));
				break;

			case TextType:
				dataSetWriter.write16Bit(nv.getValueText(), data.getMetricColumnLength(currentDataType, icol));
				break;

			case UnknownType:
				break;
			}
		}
	}

	/**
	 * Write headers.
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 * @throws UnsignedOutOfLimitsException
	 *           the unsigned out of limits exception
	 */
	private void writeHeaders() throws IOException, UnsignedOutOfLimitsException {
		writer.writeHeader();
		DataGroupWriter dataGroupWriter = writer.getDataGroupWriter(0);
		dataGroupWriter.writeHeader();

		int n = dataGroupWriter.getDataSetWriterCnt();
		for (int i = 0; i < n; i++) {
			dataSetWriter = dataGroupWriter.getDataSetWriter(i);
			dataSetWriter.writeHeader();
			entryPos.put(indexToDataType.get(i), setFilePositions());
			dataGroupWriter.updateNextDataGroupPos();
		}
	}

	/**
	 * Sets the file positions.
	 * 
	 * @return the u int
	 * 
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 * @throws UnsignedOutOfLimitsException
	 *           the unsigned out of limits exception
	 */
	private UInt setFilePositions() throws IOException, UnsignedOutOfLimitsException {
		int dataSetSz = dataSetWriter.getSize();
		UInt offset = writer.getFilePos();
		writer.seekFromCurrentPos(new UInt(0xFFFFFFFFL & (long)(dataSetSz + 1)));
		dataSetWriter.updateNextDataSetOffset();
		return offset;
	}
}
