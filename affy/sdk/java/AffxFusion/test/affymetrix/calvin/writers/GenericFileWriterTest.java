package affymetrix.calvin.writers;

import java.io.IOException;
import java.util.Arrays;

import junit.framework.TestCase;
import affymetrix.calvin.data.FileHeader;
import affymetrix.calvin.data.GenericData;
import affymetrix.calvin.data.GenericDataHeader;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parsers.GenericFileReader;
import affymetrix.calvin.parsers.InvalidFileTypeException;
import affymetrix.calvin.parsers.InvalidVersionException;
import affymetrix.calvin.utils.AffymetrixGuidType;

public class GenericFileWriterTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testWrite() throws IOException, UnsignedOutOfLimitsException, InvalidFileTypeException,
			InvalidVersionException {
		FileHeader hdr = new FileHeader();
		hdr.setFilename("generic_file_writer");
		GenericDataHeader gHdr = new GenericDataHeader();
		gHdr.setFileCreationTime("creation time");
		AffymetrixGuidType guid = new AffymetrixGuidType();
		gHdr.setFileId(guid);
		gHdr.setFileTypeId("file type id");
		gHdr.setLocale("locale");
		hdr.setGenericDataHdr(gHdr);
		GenericFileWriter writer = new GenericFileWriter(hdr);
		writer.writeHeader();
		writer.close();
		writer = null;
		hdr = null;
		gHdr = null;
		assertTrue(true);

		GenericFileReader reader = new GenericFileReader();
		reader.setFilename("generic_file_writer");
		GenericData data = new GenericData();
		reader.readHeader(data, GenericFileReader.ReadHeaderOption.ReadAllHeaders);
		hdr = data.getHeader();
		gHdr = hdr.getGenericDataHdr();
		assertTrue(Arrays.equals(gHdr.getFileId().getGuid(), guid.getGuid()));
		assertEquals(gHdr.getFileTypeId(), "file type id");
		assertEquals(gHdr.getFileCreationTime(), "creation time");
		assertEquals(gHdr.getLocale(), "locale");
	}
}
