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
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import affymetrix.calvin.data.DataGroupHeader;
import affymetrix.calvin.data.DataSetHeader;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.portability.UInt;

public class DataGroupWriter {

	BufferedFileOutputStream output;

	DataGroupHeader dataGroupHdr;

	DataGroupHeaderWriter dataGroupHdrWriter = new DataGroupHeaderWriter();

	List<DataSetWriter> dataSetWriters = new ArrayList<DataSetWriter>();

	public DataGroupWriter(BufferedFileOutputStream os, DataGroupHeader hdr) {
		output = os;
		dataGroupHdr = hdr;
		createWriters();
	}

	public DataSetWriter createDataSetWriter(DataSetHeader hdr) {
		return new DataSetWriter(output, hdr);
	}

	private void createWriters() {
		int sz = dataGroupHdr.getDataSetCnt();
		for (int i = 0; i < sz; i++) {
			DataSetWriter p = createDataSetWriter(dataGroupHdr.getDataSet(i));
			dataSetWriters.add(p);
		}
	}

	public void writeHeader() throws IOException, UnsignedOutOfLimitsException {
		dataGroupHdrWriter.write(output, dataGroupHdr);
		UInt currentPos = new UInt(output.getPosition());
		dataGroupHdrWriter.updateDataSetPos(output, currentPos);
	}

	public void close() throws IOException, UnsignedOutOfLimitsException {
		UInt currentPos = new UInt(output.getPosition());
		dataGroupHdrWriter.updateNextDataGroupPos(output, currentPos);
	}

	public Iterator<DataSetWriter> getDataSetWriterIterator() {
		return dataSetWriters.iterator();
	}

	public DataSetWriter getDataSetWriter(int index) {
		return dataSetWriters.get(index);
	}

	public int getDataSetWriterCnt() {
		return dataSetWriters.size();
	}

	public String getDataGroupName() {
		return dataGroupHdr.getName();
	}

	public void updateNextDataGroupPos() throws IOException, UnsignedOutOfLimitsException {
		UInt currentPos = new UInt(output.getPosition());
		dataGroupHdrWriter.updateNextDataGroupPos(output, currentPos);
	}
}
