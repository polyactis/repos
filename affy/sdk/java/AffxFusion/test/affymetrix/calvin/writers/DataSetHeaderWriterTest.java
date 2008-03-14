package affymetrix.calvin.writers;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import junit.framework.TestCase;
import affymetrix.calvin.data.DataSetHeader;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parsers.DataSetHeaderReader;

public class DataSetHeaderWriterTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testWrite() throws FileNotFoundException, IOException, UnsignedOutOfLimitsException {
		DataSetHeaderWriter writer = new DataSetHeaderWriter();
		String f = "data_dataGroup_header";
		BufferedFileOutputStream os = new BufferedFileOutputStream(f);
		DataSetHeader hdr = new DataSetHeader();
		hdr.setName("xyz123");
		hdr.addAsciiColumn("xyz", 128);
		hdr.addIntColumn("123");
		hdr.setRowCnt(3);
		writer.write(os, hdr);
		hdr.clear();
		os.close();
		assertTrue(true);

		// read header and verify data
		DataSetHeaderReader reader = new DataSetHeaderReader();
		FileInputStream is = new FileInputStream(f);
		reader.read(is, hdr);
		assertEquals(hdr.getName(), "xyz123");
	}
}
