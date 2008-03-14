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

import affymetrix.calvin.data.CHPQuantificationDetectionData;
import affymetrix.calvin.data.ProbeSetQuantificationDetectionData;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.portability.UInt;

public class CHPQuantificationDetectionFileWriter extends CHPFileWriterBase {

	/** The data set writer. */
	private DataSetWriter dataSetWriter;

	/** The file position of the entry for each data set. */
	private UInt entryPos;

	/** The maximum probe set name length. */
	private int maxProbeSetName;

	public CHPQuantificationDetectionFileWriter(CHPQuantificationDetectionData p) throws IOException,
			UnsignedOutOfLimitsException {
		super(p.getFileHeader());
		maxProbeSetName = p.getMaxProbeSetName();
		writeHeaders();
	}

	public void writeHeaders() throws IOException, UnsignedOutOfLimitsException {
		writer.writeHeader();
		DataGroupWriter dataGroupWriter = writer.getDataGroupWriter(0);
		dataGroupWriter.writeHeader();
		Iterator<DataSetWriter> it = dataGroupWriter.getDataSetWriterIterator();
		if (it.hasNext()) {
			dataSetWriter = it.next();
			dataSetWriter.writeHeader();
			entryPos = setFilePositions();
			dataGroupWriter.updateNextDataGroupPos();
		}
	}

	public void seekToDataSet() throws IOException {
		DataGroupWriter dataGroupWriter = writer.getDataGroupWriter(0);
		dataSetWriter = dataGroupWriter.getDataSetWriter(0);
		writer.seekFromBeginPos(entryPos);
	}

	public void writeEntry(ProbeSetQuantificationDetectionData p) throws IOException {
		if (maxProbeSetName == -1) {
			dataSetWriter.write(p.getId());
		}
		else {
			dataSetWriter.write8Bit(p.getName(), maxProbeSetName);
		}
		dataSetWriter.write(p.getQuantification());
		dataSetWriter.write(p.getPValue());
	}

	public UInt setFilePositions() throws IOException, UnsignedOutOfLimitsException {
		int dataSetSz = dataSetWriter.getSize();
		UInt offset = writer.getFilePos();
		writer.seekFromCurrentPos(new UInt(0xFFFFFFFFL & (long)(dataSetSz + 1)));
		dataSetWriter.updateNextDataSetOffset();
		return offset;
	}
}
