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

import affymetrix.calvin.data.ColumnInfo;
import affymetrix.calvin.data.DataSetHeader;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parameter.ParameterNameValue;
import affymetrix.portability.UByte;
import affymetrix.portability.UInt;

public class DataSetHeaderWriter {

	private UInt dataPos;

	private UInt nextDataSetPos;

	public DataSetHeaderWriter() {
		dataPos = new UInt();
		nextDataSetPos = new UInt();
	}

	public void UpdateDataOffset(BufferedFileOutputStream output, UInt pos) throws IOException {
		if (dataPos.greaterThan(0)) {
			output.setPosition(dataPos.toLong());
			FileOutput.writeUInt32(output, pos);
			output.setPosition(pos.toLong());
			// output.flush();
		}
	}

	public void updateNextDataSetOffset(BufferedFileOutputStream output, UInt pos) throws IOException {
		if (nextDataSetPos.greaterThan(0)) {
			output.setPosition(nextDataSetPos.toLong());
			FileOutput.writeUInt32(output, pos);
			output.setPosition(pos.toLong());
			// output.flush();
		}
	}

	public void write(BufferedFileOutputStream output, DataSetHeader hdr) throws IOException,
			UnsignedOutOfLimitsException {
		writeDataOffset(output, UInt.ZERO);
		writeNextDataSetOffset(output, UInt.ZERO);
		writeName(output, hdr);
		writeNameValCnt(output, hdr);
		writeNameValParams(output, hdr);
		writeColumnCnt(output, hdr);
		writeColumnTypes(output, hdr);
		writeRowCnt(output, hdr);
	}

	private void writeColumnCnt(BufferedFileOutputStream output, DataSetHeader hdr) throws IOException {
		FileOutput.writeInt32(output, hdr.getColumnCnt());
	}

	private void writeColumnTypes(BufferedFileOutputStream output, DataSetHeader hdr) throws IOException,
			UnsignedOutOfLimitsException {
		// write the types
		for (int i = 0; i < hdr.getColumnCnt(); i++) {
			ColumnInfo col = hdr.getColumnInfo(i);
			FileOutput.writeString16(output, col.getName());
			FileOutput.writeUInt8(output, new UByte((short)col.getColumnType().ordinal()));
			FileOutput.writeInt32(output, col.getSize());
		}
	}

	private void writeDataOffset(BufferedFileOutputStream output, UInt pos) throws IOException,
			UnsignedOutOfLimitsException {
		dataPos.set(output.getPosition());
		FileOutput.writeUInt32(output, pos);
	}

	private void writeName(BufferedFileOutputStream output, DataSetHeader hdr) throws IOException {
		FileOutput.writeString16(output, hdr.getName());
	}

	private void writeNameValCnt(BufferedFileOutputStream output, DataSetHeader hdr) throws IOException {
		FileOutput.writeInt32(output, hdr.getNameValParamCnt());
	}

	private void writeNameValParams(BufferedFileOutputStream output, DataSetHeader hdr) throws IOException {
		Iterator i = hdr.getNameValParameters().iterator();
		while (i.hasNext()) {
			ParameterNameValue p = (ParameterNameValue)i.next();
			FileOutput.writeString16(output, p.getName());
			FileOutput.writeBlob(output, p.getMIMEValue().getBytes());
			FileOutput.writeString16(output, p.getMIMEType());
		}
	}

	private void writeNextDataSetOffset(BufferedFileOutputStream output, UInt pos) throws IOException,
			UnsignedOutOfLimitsException {
		nextDataSetPos.set(output.getPosition());
		FileOutput.writeUInt32(output, pos);
	}

	private void writeRowCnt(BufferedFileOutputStream output, DataSetHeader hdr) throws IOException {
		FileOutput.writeInt32(output, hdr.getRowCnt());
	}
}
