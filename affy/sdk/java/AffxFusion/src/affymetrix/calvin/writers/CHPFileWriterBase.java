package affymetrix.calvin.writers;

import java.io.FileNotFoundException;
import java.io.IOException;

import affymetrix.calvin.data.FileHeader;

public class CHPFileWriterBase {

	/** The file writer. */
	protected GenericFileWriter writer;

	CHPFileWriterBase(FileHeader hdr) throws IOException, FileNotFoundException {
		writer = new GenericFileWriter(hdr);
	}

	public void flush() throws IOException {
		writer.flush();
	}

	public void close() throws IOException {
		writer.close();
	}
}
