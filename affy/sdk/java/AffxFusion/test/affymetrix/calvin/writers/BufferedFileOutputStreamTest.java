package affymetrix.calvin.writers;

import java.io.File;
import java.io.FileOutputStream;

import junit.framework.TestCase;
import affymetrix.calvin.utils.IOUtils;
import affymetrix.portability.UInt;

public class BufferedFileOutputStreamTest extends TestCase {

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testWrite() throws Exception {
		byte[] msg = "test message for buffered output".getBytes(IOUtils.ASCII_CHARSET);
		byte[] num = new UInt(362436).getBytes();
		File bufFile = new File("BufferedFileOuput");
		BufferedFileOutputStream buf = new BufferedFileOutputStream(bufFile);
		long begin = System.currentTimeMillis();
		for (int i = 0; i < 100000; i++) {
			buf.write(msg);
		}
		long bufPos1 = buf.getPosition();
		buf.setPosition(bufPos1 + 100);
		for (int i = 0; i < 100000; i++) {
			buf.write(num);
		}
		long bufPos2 = buf.getPosition();
		buf.setPosition(bufPos2 - 500);
		for (int i = 0; i < 100000; i++) {
			buf.write(num);
		}
		buf.flush();
		long end = System.currentTimeMillis();
		long bufElapsed = end - begin;

		File osFile = new File("UnbufferedFileOuput");
		FileOutputStream os = new FileOutputStream(osFile);
		begin = System.currentTimeMillis();
		for (int i = 0; i < 100000; i++) {
			os.write(msg);
		}
		long osPos1 = os.getChannel().position();
		os.getChannel().position(osPos1 + 100);
		for (int i = 0; i < 100000; i++) {
			os.write(num);
		}
		long osPos2 = os.getChannel().position();
		os.getChannel().position(osPos2 - 500);
		for (int i = 0; i < 100000; i++) {
			os.write(num);
		}
		os.flush();
		end = System.currentTimeMillis();
		long osElapsed = end - begin;
		System.out.println("Buffered: " + bufElapsed + ", Unbuffered: " + osElapsed);
		assertTrue(osElapsed > bufElapsed);
		assertEquals(os.getChannel().position(), buf.getPosition());
		assertEquals(osPos1, bufPos1);
		assertEquals(osPos2, bufPos2);

		bufFile.delete();
		osFile.delete();
	}
}
