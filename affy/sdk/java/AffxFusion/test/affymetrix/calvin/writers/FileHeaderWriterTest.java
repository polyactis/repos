package affymetrix.calvin.writers;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import junit.framework.TestCase;
import affymetrix.calvin.data.FileHeader;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parsers.FileHeaderReader;
import affymetrix.calvin.parsers.InvalidFileTypeException;
import affymetrix.calvin.parsers.InvalidVersionException;

public class FileHeaderWriterTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testWrite() throws FileNotFoundException, IOException, UnsignedOutOfLimitsException,
			InvalidVersionException, InvalidFileTypeException {
		FileHeaderWriter writer = new FileHeaderWriter();
		String f = "file_header";
		BufferedFileOutputStream os = new BufferedFileOutputStream(f, true);
		FileHeader fHdr = new FileHeader();
		fHdr.setFilename(f);
		writer.write(os, fHdr);
		os.close();
		fHdr.clear();
		assertTrue(true);

		// read header and verify data
		FileInputStream is = new FileInputStream(f);
		FileHeaderReader reader = new FileHeaderReader(is, fHdr);
		reader.read();
		assertEquals(fHdr.getMagicNumber(), FileHeader.MAGIC_NUM);
		assertEquals(fHdr.getVersion(), FileHeader.VERSION);
	}
}
