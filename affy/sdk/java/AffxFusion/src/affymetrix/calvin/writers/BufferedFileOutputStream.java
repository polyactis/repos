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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;

public class BufferedFileOutputStream {

	private static final int MAX_BUFFER_SIZE = 5242880; // 5 MB

	FileOutputStream os = null;

	byte[] buf = null;

	int current = 0;

	int maxSz = 0;

	public BufferedFileOutputStream(String filename) throws FileNotFoundException {
		this(filename, false);
	}

	public BufferedFileOutputStream(String filename, boolean append) throws FileNotFoundException {
		this(filename, append, MAX_BUFFER_SIZE);
	}

	public BufferedFileOutputStream(String filename, boolean append, int size) throws FileNotFoundException {
		this(new File(filename), append, size);
	}

	public BufferedFileOutputStream(File file) throws FileNotFoundException {
		this(file, false);
	}

	public BufferedFileOutputStream(File file, boolean append) throws FileNotFoundException {
		this(file, append, MAX_BUFFER_SIZE);
	}

	public BufferedFileOutputStream(File file, boolean append, int size) throws FileNotFoundException {
		this(new FileOutputStream(file, append), size);
	}

	public BufferedFileOutputStream(FileOutputStream fos) throws FileNotFoundException {
		this(fos, MAX_BUFFER_SIZE);
	}

	public BufferedFileOutputStream(FileOutputStream fos, int size) throws FileNotFoundException {
		os = fos;
		buf = new byte[size];
		maxSz = size;
	}

	public void flush() throws IOException {
		if ((maxSz > 0) && (current > 0)) {
			os.write(buf, 0, current);
			current = 0;
		}
	}

	public void close() throws IOException {
		flush();
		os.close();
	}

	public long getSize() throws IOException {
		return os.getChannel().size();
	}

	public long getPosition() throws IOException {
		return os.getChannel().position() + current;
	}

	public void setPosition(long p) throws IOException {
		flush();
		os.getChannel().position(p);
	}

	public void write(byte b) throws IOException {
		if (maxSz > 0) {
			buffer(b);
		}
		else {
			os.write(b);
		}
	}

	public void write(byte[] b) throws IOException {
		write(b, 0, b.length);
	}

	public void write(byte[] b, int offset, int length) throws IOException {
		if (maxSz > 0) {
			for (int i = offset; i < length; i++) {
				buffer(b[i]);
			}
		}
		else {
			os.write(b, offset, length);
		}
	}

	private void buffer(byte b) throws IOException {
		if (current >= maxSz) {
			flush();
		}
		buf[current] = b;
		current++;
	}

	@Override
	protected void finalize() throws Throwable {
		close();
		super.finalize();
	}
}
