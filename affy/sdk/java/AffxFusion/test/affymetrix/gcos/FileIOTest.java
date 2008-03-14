/*
 * FileIOTest.java
 * JUnit based test
 *
 * Created on October 13, 2005, 8:02 AM
 */

package affymetrix.gcos;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * 
 * @author ljevon
 */
public class FileIOTest extends TestCase {

	public FileIOTest(String testName) {
		super(testName);
	}

	protected void setUp() throws Exception {
	}

	protected void tearDown() throws Exception {
	}

	public static Test suite() {
		TestSuite suite = new TestSuite(FileIOTest.class);

		return suite;
	}

	public void testReadNextLine() throws Exception {
		String prefix = new File("..\\..\\..\\CPPTest\\data\\test.file").getCanonicalPath();
		FileReader f = null;
		BufferedReader b = null;
		try {
			f = new FileReader(prefix + ".unix");
			b = new BufferedReader(f);
		}
		catch (Throwable t) {
		}
		String str = FileIO.readNextLine(b);
		assertEquals(str, "line 1");
		str = FileIO.readNextLine(b);
		assertEquals(str, "line 2");
		str = FileIO.readNextLine(b);
		assertEquals(str, "line 3");
		str = FileIO.readNextLine(b);
		assertEquals(str, null);

		try {
			f = new FileReader(prefix + ".dos");
			b = new BufferedReader(f);
		}
		catch (Throwable t) {
		}
		str = FileIO.readNextLine(b);
		assertEquals(str, "line 1");
		str = FileIO.readNextLine(b);
		assertEquals(str, "line 2");
		str = FileIO.readNextLine(b);
		assertEquals(str, "line 3");
		str = FileIO.readNextLine(b);
		assertEquals(str, null);

	}

	public void testReadInt_T() throws Exception {
		String prefix = new File("..\\..\\..\\CPPTest\\data\\test.file").getCanonicalPath();
		FileInputStream fis = null;
		try {
			fis = new FileInputStream(prefix + ".int");
		}
		catch (Throwable t) {
		}
		int ival = FileIO.readInt32(fis);
		assertEquals(ival, 10);
		ival = FileIO.readInt32(fis);
		assertEquals(ival, 20);
		ival = FileIO.readInt32(fis);
		assertEquals(ival, 30);

		try {
			fis = new FileInputStream(prefix + ".uint");
		}
		catch (Throwable t) {
		}
		ival = FileIO.readInt32(fis);
		assertEquals(ival, 101);
		ival = FileIO.readInt32(fis);
		assertEquals(ival, 202);
		ival = FileIO.readInt32(fis);
		assertEquals(ival, 303);
	}

	public void testReadFloat() throws Exception {
		String prefix = new File("..\\..\\..\\CPPTest\\data\\test.file").getCanonicalPath();
		FileInputStream fis = null;
		try {
			fis = new FileInputStream(prefix + ".float");
		}
		catch (Throwable t) {
		}
		float fval = FileIO.readFloat(fis);
		assertEquals(fval, 321.123f, 0.000001f);
		fval = FileIO.readFloat(fis);
		assertEquals(fval, 123.123f, 0.000001f);
		fval = FileIO.readFloat(fis);
		assertEquals(fval, 543.56f, 0.000001f);
	}

	public void testReadString() throws Exception {
		String prefix = new File("..\\..\\..\\CPPTest\\data\\test.file").getCanonicalPath();
		FileInputStream fis = null;
		try {
			fis = new FileInputStream(prefix + ".string");
		}
		catch (Throwable t) {
		}
		String sval = FileIO.readString(fis);
		assertEquals(sval, "string");
		sval = FileIO.readString(fis);
		assertEquals(sval, "test");
		sval = FileIO.readString(fis);
		assertEquals(sval, "case");
	}
}
