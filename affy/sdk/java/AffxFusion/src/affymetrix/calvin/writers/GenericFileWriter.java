//////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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
///////////////////////////////////////////////////////////////////////////////////
package affymetrix.calvin.writers;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import affymetrix.calvin.data.DataGroupHeader;
import affymetrix.calvin.data.FileHeader;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.portability.UInt;

public class GenericFileWriter {

	private BufferedFileOutputStream output = null;

	private FileHeader fileHdr = null;

	private List<DataGroupWriter> dataGrpWriters = new ArrayList<DataGroupWriter>();

	public GenericFileWriter(FileHeader hdr) throws FileNotFoundException, IOException {
		this(hdr, true);
	}

	/**
	 * Instantiates a new generic file writer.
	 * 
	 * @param hdr
	 *          file header
	 * @param truncate
	 *          truncate file if true
	 * 
	 * @throws FileNotFoundException
	 *           file not found exception
	 * @throws IOException
	 *           Signals that an I/O exception has occurred.
	 */
	public GenericFileWriter(FileHeader hdr, boolean truncate) throws FileNotFoundException, IOException {
		fileHdr = hdr;
		openFileOStream(fileHdr.getFilename(), truncate);
		createWriters();
	}

	@Override
	protected void finalize() throws Throwable {
		if (dataGrpWriters != null) {
			dataGrpWriters.clear();
		}
		if (output != null) {
			output.flush();
			output.close();
		}
		super.finalize();
	}

	public void flush() throws IOException {
		if (output != null) {
			output.flush();
		}
	}

	public void close() throws IOException {
		if (output != null) {
			output.close();
		}
	}

	public UInt getFilePos() throws IOException, UnsignedOutOfLimitsException {
		return new UInt(output.getPosition());
	}

	private DataGroupWriter createDataGroupWriter(DataGroupHeader hdr) {
		return new DataGroupWriter(output, hdr);
	}

	private void createWriters() {
		int sz = fileHdr.getDataGroupCnt();
		for (int i = 0; i < sz; i++) {
			DataGroupWriter p = createDataGroupWriter(fileHdr.getDataGroup(i));
			dataGrpWriters.add(p);
		}
	}

	public DataGroupWriter getDataGroupWriter(int index) {
		return dataGrpWriters.get(index);
	}

	public Iterator<DataGroupWriter> getDataGroupWriterIterator() {
		return dataGrpWriters.iterator();
	}

	public void openFileOStream(String file, boolean truncate) throws FileNotFoundException, IOException {
		if (truncate) {
			output = new BufferedFileOutputStream(file);
		}
		else {
			output = new BufferedFileOutputStream(file, true);
		}
	}

	public void seekFromBeginPos(UInt offset) throws IOException {
		output.setPosition(offset.toLong());
	}

	public void seekFromCurrentPos(UInt offset) throws IOException, UnsignedOutOfLimitsException {
		long pos = output.getPosition();
		output.setPosition(offset.add(pos));
	}

	public void seekFromEndPos(UInt offset) throws IOException {
		long pos = output.getSize() - offset.toLong();
		output.setPosition(pos);
	}

	public void write(UInt p) throws IOException {
		FileOutput.writeUInt32(output, p);
	}

	public void writeHeader() throws IOException, UnsignedOutOfLimitsException {
		FileHeaderWriter fileHdrWriter = new FileHeaderWriter();
		fileHdrWriter.write(output, fileHdr);
		fileHdrWriter.updateDataGroupOffset(output, new UInt(output.getPosition()));
	}
}
