package affymetrix.calvin.writers;

import junit.framework.TestCase;
import affymetrix.calvin.data.GenericDataHeader;
import affymetrix.calvin.utils.AffymetrixGuidType;

public class GenericDataHeaderWriterTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testWrite() throws Exception {
		String f = "generic_data_header";
		BufferedFileOutputStream os = new BufferedFileOutputStream(f);
		GenericDataHeaderWriter writer = new GenericDataHeaderWriter();
		GenericDataHeader g = new GenericDataHeader();
		AffymetrixGuidType guid = new AffymetrixGuidType();
		g.setFileId(guid);
		writer.write(os, g);
		os.close();
		assertTrue(true);

		// TODO: use generic data header reader to read header and verify data
	}
}
