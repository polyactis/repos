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

import affymetrix.calvin.data.DataSetHeader;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.portability.UByte;
import affymetrix.portability.UInt;
import affymetrix.portability.UShort;

public class DataSetWriter {

	private DataSetHeader dataSetHdr;

	private DataSetHeaderWriter dataSetHdrWriter;

	private BufferedFileOutputStream output;

	public DataSetWriter(BufferedFileOutputStream bw, DataSetHeader hdr) {
		output = bw;
		dataSetHdr = hdr;
		dataSetHdrWriter = new DataSetHeaderWriter();
	}

	public String getName() {
		return dataSetHdr.getName();
	}

	public int getSize() {
		return dataSetHdr.getDataSize();
	}

	public void writeHeader() throws IOException, UnsignedOutOfLimitsException {
		dataSetHdrWriter.write(output, dataSetHdr);
		UInt currentPos = new UInt(output.getPosition());
		dataSetHdr.setDataStartFilePos(currentPos);
		dataSetHdrWriter.UpdateDataOffset(output, currentPos);
	}

	public void updateNextDataSetOffset() throws IOException, UnsignedOutOfLimitsException {
		UInt currentPos = new UInt(output.getPosition());
		dataSetHdr.setNextSetFilePos(currentPos);
		dataSetHdrWriter.updateNextDataSetOffset(output, currentPos);
	}

	public void write(byte p) throws IOException {
		FileOutput.writeInt8(output, p);
	}

	public void write(UByte p) throws IOException {
		FileOutput.writeUInt8(output, p);
	}

	public void write(short p) throws IOException {
		FileOutput.writeInt16(output, p);
	}

	public void write(UShort p) throws IOException {
		FileOutput.writeUInt16(output, p);
	}

	public void write(int p) throws IOException {
		FileOutput.writeInt32(output, p);
	}

	public void write(UInt p) throws IOException {
		FileOutput.writeUInt32(output, p);
	}

	public void write(float p) throws IOException {
		FileOutput.writeFloat(output, p);
	}

	public void write8Bit(String p, int maxLn) throws IOException {
		FileOutput.writeString8(output, p, maxLn);
	}

	public void write16Bit(String p, int maxLn) throws IOException {
		FileOutput.writeString16(output, p, maxLn);
	}

}
