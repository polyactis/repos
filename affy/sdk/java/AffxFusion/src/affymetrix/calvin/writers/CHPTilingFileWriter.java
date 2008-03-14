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

import affymetrix.calvin.data.CHPTilingData;
import affymetrix.calvin.data.CHPTilingEntry;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.portability.UInt;

public class CHPTilingFileWriter extends CHPFileWriterBase {

	/** The data set writer. */
	private DataSetWriter dataSetWriter;

	/** The file position of the entry for each data set. */
	private UInt[] entryPos;

	CHPTilingFileWriter(CHPTilingData p) throws IOException, UnsignedOutOfLimitsException {
		super(p.getFileHeader());
		entryPos = new UInt[p.getNumberSequences()];
		writeHeaders();
	}

	public void seekToDataSet(int index) throws IOException {
		DataGroupWriter dataGroupWriter = writer.getDataGroupWriter(0);
		dataSetWriter = dataGroupWriter.getDataSetWriter(index);
		writer.seekFromBeginPos(entryPos[index]);
		// dataSetIndex = index;
	}

	public void writeTilingEntry(CHPTilingEntry p) throws IOException {
		dataSetWriter.write(p.getPosition());
		dataSetWriter.write(p.getValue());
	}

	private UInt setFilePositions() throws IOException, UnsignedOutOfLimitsException {
		int dataSetSz = dataSetWriter.getSize();
		UInt offset = writer.getFilePos();
		writer.seekFromCurrentPos(new UInt(0xFFFFFFFFL & (long)(dataSetSz + 1)));
		dataSetWriter.updateNextDataSetOffset();
		return offset;
	}

	private void writeHeaders() throws IOException, UnsignedOutOfLimitsException {
		writer.writeHeader();
		DataGroupWriter dataGroupWriter = writer.getDataGroupWriter(0);
		dataGroupWriter.writeHeader();

		int iSet = 0;
		Iterator<DataSetWriter> it = dataGroupWriter.getDataSetWriterIterator();
		while (it.hasNext()) {
			dataSetWriter = it.next();
			dataSetWriter.writeHeader();
			entryPos[iSet++] = setFilePositions();
		}
		dataGroupWriter.updateNextDataGroupPos();
	}
}
