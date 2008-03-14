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

import affymetrix.calvin.data.DataGroupHeader;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.portability.UInt;

public class DataGroupHeaderWriter {

	private UInt dataSetPos = new UInt();

	private UInt nextDataGroupPos = new UInt();

	public DataGroupHeaderWriter() {
	}

	public void write(BufferedFileOutputStream output, DataGroupHeader hdr) throws IOException,
			UnsignedOutOfLimitsException {
		writeNextDataGroupPos(output, UInt.ZERO);
		writeDataSetPos(output, UInt.ZERO);
		writeDataSetCnt(output, hdr);
		writeName(output, hdr);
		output.flush();
	}

	public void updateNextDataGroupPos(BufferedFileOutputStream output, UInt pos) throws IOException {
		if (nextDataGroupPos.greaterThan(0)) {
			output.setPosition(nextDataGroupPos.toLong());
			FileOutput.writeUInt32(output, pos);
			output.setPosition(pos.toLong());
			// output.flush();
		}
	}

	public void updateDataSetPos(BufferedFileOutputStream output, UInt pos) throws IOException {
		if (dataSetPos.greaterThan(0)) {
			output.setPosition(dataSetPos.toLong());
			FileOutput.writeUInt32(output, pos);
			output.setPosition(pos.toLong());
			// output.flush();
		}
	}

	private void writeDataSetPos(BufferedFileOutputStream output, UInt pos) throws IOException,
			UnsignedOutOfLimitsException {
		dataSetPos.set(output.getPosition());
		FileOutput.writeUInt32(output, pos);
	}

	private void writeNextDataGroupPos(BufferedFileOutputStream output, UInt pos) throws IOException,
			UnsignedOutOfLimitsException {
		nextDataGroupPos.set(output.getPosition());
		FileOutput.writeUInt32(output, pos);
	}

	private void writeName(BufferedFileOutputStream output, DataGroupHeader hdr) throws IOException {
		FileOutput.writeString16(output, hdr.getName());
	}

	private void writeDataSetCnt(BufferedFileOutputStream output, DataGroupHeader hdr) throws IOException {
		FileOutput.writeInt32(output, hdr.getDataSetCnt());
	}

}
