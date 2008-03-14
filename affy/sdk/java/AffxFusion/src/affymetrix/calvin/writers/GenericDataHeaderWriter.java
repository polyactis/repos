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

import affymetrix.calvin.data.GenericDataHeader;
import affymetrix.calvin.parameter.ParameterNameValue;

public class GenericDataHeaderWriter {

	public GenericDataHeaderWriter() {
	}

	public void write(BufferedFileOutputStream os, GenericDataHeader g) throws IOException {
		writeFileTypeId(os, g);
		writeFileId(os, g);
		writeFileCreationTime(os, g);
		writeLocale(os, g);
		writeNameValParamCnt(os, g);
		writeNameValParams(os, g);
		writeParentHdrCnt(os, g);
		writeParentHdrs(os, g);
	}

	private void writeFileTypeId(BufferedFileOutputStream os, GenericDataHeader g) throws IOException {
		FileOutput.writeString8(os, g.getFileTypeId());
	}

	private void writeFileId(BufferedFileOutputStream os, GenericDataHeader g) throws IOException {
		FileOutput.writeBlob(os, g.getFileId().getGuid());
	}

	private void writeFileCreationTime(BufferedFileOutputStream os, GenericDataHeader g) throws IOException {
		FileOutput.writeString16(os, g.getFileCreationTime());
	}

	private void writeLocale(BufferedFileOutputStream os, GenericDataHeader g) throws IOException {
		FileOutput.writeString16(os, g.getLocale());
	}

	private void writeNameValParamCnt(BufferedFileOutputStream os, GenericDataHeader g) throws IOException {
		FileOutput.writeInt32(os, g.getNameValParamCnt());
	}

	private void writeNameValParams(BufferedFileOutputStream os, GenericDataHeader g) throws IOException {
		Iterator<ParameterNameValue> i = g.getNameValParams().iterator();
		while (i.hasNext()) {
			ParameterNameValue p = i.next();
			FileOutput.writeString16(os, p.getName());
			FileOutput.writeBlob(os, p.getMIMEValue().getBytes());
			FileOutput.writeString16(os, p.getMIMEType());
		}
	}

	private void writeParentHdrCnt(BufferedFileOutputStream os, GenericDataHeader g) throws IOException {
		FileOutput.writeInt32(os, g.getParentCnt());
	}

	private void writeParentHdrs(BufferedFileOutputStream os, GenericDataHeader g) throws IOException {
		Iterator<GenericDataHeader> i = g.getParents().iterator();
		while (i.hasNext()) {
			write(os, i.next());
		}
	}
}
