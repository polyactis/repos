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

import affymetrix.calvin.data.FileHeader;
import affymetrix.calvin.data.GenericDataHeader;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.utils.AffymetrixGuidType;
import affymetrix.portability.UInt;

public class FileHeaderWriter {

	private UInt dataSetOffsetPos = new UInt();

	public FileHeaderWriter() {
	}

	public void write(BufferedFileOutputStream os, FileHeader hdr) throws IOException, UnsignedOutOfLimitsException {
		writeMagicNumber(os, hdr);
		writeVersion(os, hdr);
		writeDataGroupCnt(os, hdr);
		writeDataGroupOffset(os, UInt.ZERO);
		writeGenericDataHdr(os, hdr);
	}

	private void writeMagicNumber(BufferedFileOutputStream os, FileHeader hdr) throws IOException {
		FileOutput.writeInt8(os, hdr.getMagicNumber());
	}

	private void writeVersion(BufferedFileOutputStream os, FileHeader g) throws IOException {
		FileOutput.writeInt8(os, g.getVersion());
	}

	private void writeDataGroupCnt(BufferedFileOutputStream os, FileHeader g) throws IOException {
		FileOutput.writeInt32(os, g.getDataGroupCnt());
	}

	private void writeDataGroupOffset(BufferedFileOutputStream os, UInt offset) throws IOException,
			UnsignedOutOfLimitsException {
		dataSetOffsetPos.set(os.getPosition());
		FileOutput.writeUInt32(os, offset);
	}

	public void updateDataGroupOffset(BufferedFileOutputStream os, UInt offset) throws IOException {
		if (dataSetOffsetPos.greaterThan(0)) {
			os.setPosition(dataSetOffsetPos.toLong());
			FileOutput.writeUInt32(os, offset);
			os.setPosition(offset.toLong());
		}
	}

	void writeGenericDataHdr(BufferedFileOutputStream os, FileHeader g) throws IOException {
		GenericDataHeaderWriter gdhWriter = new GenericDataHeaderWriter();

		// Check if a file ID has been assign, if not assign one.
		GenericDataHeader gdh = g.getGenericDataHdr();
		AffymetrixGuidType id = gdh.getFileId();
		if (id.isEmpty()) {
			id.generateGuid();
			gdh.setFileId(id);
		}
		gdhWriter.write(os, gdh);
	}

}
