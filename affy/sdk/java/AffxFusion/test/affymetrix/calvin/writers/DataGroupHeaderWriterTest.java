package affymetrix.calvin.writers;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import junit.framework.TestCase;
import affymetrix.calvin.data.DataGroupHeader;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parsers.DataGroupHeaderReader;

public class DataGroupHeaderWriterTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testWrite() throws FileNotFoundException, IOException, UnsignedOutOfLimitsException {
		DataGroupHeaderWriter writer = new DataGroupHeaderWriter();
		String f = "DataGroupHeader";
		BufferedFileOutputStream os = new BufferedFileOutputStream(f);
		DataGroupHeader grpHdr = new DataGroupHeader("DataGroup");
		writer.write(os, grpHdr);
		grpHdr.clear();
		os.close();
		assertTrue(true);

		// read header and verify data
		DataGroupHeaderReader reader = new DataGroupHeaderReader();
		FileInputStream is = new FileInputStream(f);
		reader.read(is, grpHdr);
		assertEquals(grpHdr.getName(), "DataGroup");
	}

}
