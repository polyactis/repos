package affymetrix.calvin.writers;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import junit.framework.TestCase;
import affymetrix.calvin.exception.UnsignedOutOfLimitsException;
import affymetrix.calvin.parsers.FileInput;
import affymetrix.portability.UByte;
import affymetrix.portability.UInt;
import affymetrix.portability.UShort;

public class FileOutputTest extends TestCase {

	private static final String MP = "Michelle Pfeiffer";

	private static final String JL = "Jessica Lange";

	private static final String KP = "Kelly Preston";

	protected void setUp() throws Exception {
		super.setUp();
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	public void testWriteString8() throws FileNotFoundException, IOException {
		File f = new File("./TestFile.txt");
		BufferedFileOutputStream output = new BufferedFileOutputStream(f);
		FileOutput.writeString8(output, MP);
		FileOutput.writeString8(output, JL);
		FileOutput.writeString8(output, KP);
		output.close();

		FileInputStream input = new FileInputStream(f);
		String s1 = FileInput.readString8(input);
		assertEquals(s1, MP);
		String s2 = FileInput.readString8(input);
		assertEquals(s2, JL);
		String s3 = FileInput.readString8(input);
		assertEquals(s3, KP);
		input.close();
		f.delete();
	}

	public void testWriteString16() throws FileNotFoundException, IOException {
		File f = new File("./TestFile.txt");
		BufferedFileOutputStream output = new BufferedFileOutputStream(f);
		FileOutput.writeString16(output, MP);
		FileOutput.writeString16(output, JL);
		FileOutput.writeString16(output, KP);
		output.close();

		FileInputStream input = new FileInputStream(f);
		String s1 = FileInput.readString16(input);
		assertEquals(s1, MP);
		String s2 = FileInput.readString16(input);
		assertEquals(s2, JL);
		String s3 = FileInput.readString16(input);
		assertEquals(s3, KP);
		input.close();
		f.delete();
	}

	public void testWriteBlob() throws FileNotFoundException, IOException {
		File f = new File("./TestFile.txt");
		BufferedFileOutputStream output = new BufferedFileOutputStream(f);
		FileOutput.writeBlob(output, MP.getBytes("US-ASCII"));
		FileOutput.writeBlob(output, JL.getBytes("US-ASCII"));
		FileOutput.writeBlob(output, KP.getBytes("US-ASCII"));
		output.close();

		FileInputStream input = new FileInputStream(f);
		byte[] b1 = FileInput.readBlob(input);
		assertEquals(new String(b1, "US-ASCII"), MP);
		byte[] b2 = FileInput.readBlob(input);
		assertEquals(new String(b2, "US-ASCII"), JL);
		byte[] b3 = FileInput.readBlob(input);
		assertEquals(new String(b3, "US-ASCII"), KP);
		input.close();
		f.delete();
	}

	public void testWriteUInt32() throws FileNotFoundException, IOException, UnsignedOutOfLimitsException {
		File f = new File("./TestFile.txt");
		BufferedFileOutputStream output = new BufferedFileOutputStream(f);
		UInt i1 = new UInt(2334567);
		FileOutput.writeUInt32(output, i1);
		UInt i2 = new UInt(577896);
		FileOutput.writeUInt32(output, i2);
		UInt i3 = new UInt(3789654);
		FileOutput.writeUInt32(output, i3);
		output.close();

		FileInputStream input = new FileInputStream(f);
		UInt b1 = FileInput.readUInt32(input);
		assertEquals(b1.toLong(), i1.toLong());
		UInt b2 = FileInput.readUInt32(input);
		assertEquals(b2.toLong(), i2.toLong());
		UInt b3 = FileInput.readUInt32(input);
		assertEquals(b3.toLong(), i3.toLong());
		input.close();
		f.delete();
	}

	public void testWriteInt32() throws FileNotFoundException, IOException, UnsignedOutOfLimitsException {
		File f = new File("./TestFile.txt");
		BufferedFileOutputStream output = new BufferedFileOutputStream(f);
		int i1 = 2334567;
		FileOutput.writeInt32(output, i1);
		int i2 = 577896;
		FileOutput.writeInt32(output, i2);
		int i3 = 3789654;
		FileOutput.writeInt32(output, i3);
		output.close();

		FileInputStream input = new FileInputStream(f);
		int b1 = FileInput.readInt32(input);
		assertEquals(b1, i1);
		int b2 = FileInput.readInt32(input);
		assertEquals(b2, i2);
		int b3 = FileInput.readInt32(input);
		assertEquals(b3, i3);
		input.close();
		f.delete();
	}

	public void testWriteInt16() throws FileNotFoundException, IOException, UnsignedOutOfLimitsException {
		File f = new File("./TestFile.txt");
		BufferedFileOutputStream output = new BufferedFileOutputStream(f);
		short i1 = (short)4567;
		FileOutput.writeInt16(output, i1);
		short i2 = (short)7896;
		FileOutput.writeInt16(output, i2);
		short i3 = (short)9654;
		FileOutput.writeInt16(output, i3);
		output.close();

		FileInputStream input = new FileInputStream(f);
		short b1 = FileInput.readInt16(input);
		assertEquals(b1, i1);
		short b2 = FileInput.readInt16(input);
		assertEquals(b2, i2);
		short b3 = FileInput.readInt16(input);
		assertEquals(b3, i3);
		input.close();
		f.delete();
	}

	public void testWriteUInt16() throws FileNotFoundException, IOException, UnsignedOutOfLimitsException {
		File f = new File("./TestFile.txt");
		BufferedFileOutputStream output = new BufferedFileOutputStream(f);
		UShort i1 = new UShort(34567);
		FileOutput.writeUInt16(output, i1);
		UShort i2 = new UShort(57896);
		FileOutput.writeUInt16(output, i2);
		UShort i3 = new UShort(9654);
		FileOutput.writeUInt16(output, i3);
		output.close();

		FileInputStream input = new FileInputStream(f);
		UShort b1 = FileInput.readUInt16(input);
		assertEquals(b1.toInt(), i1.toInt());
		UShort b2 = FileInput.readUInt16(input);
		assertEquals(b2.toInt(), i2.toInt());
		UShort b3 = FileInput.readUInt16(input);
		assertEquals(b3.toInt(), i3.toInt());
		input.close();
		f.delete();
	}

	public void testWriteInt8() throws FileNotFoundException, IOException, UnsignedOutOfLimitsException {
		File f = new File("./TestFile.txt");
		BufferedFileOutputStream output = new BufferedFileOutputStream(f);
		byte i1 = (byte)67;
		FileOutput.writeInt8(output, i1);
		byte i2 = (byte)96;
		FileOutput.writeInt8(output, i2);
		byte i3 = (byte)54;
		FileOutput.writeInt8(output, i3);
		output.close();

		FileInputStream input = new FileInputStream(f);
		byte b1 = FileInput.readInt8(input);
		assertEquals(b1, i1);
		byte b2 = FileInput.readInt8(input);
		assertEquals(b2, i2);
		byte b3 = FileInput.readInt8(input);
		assertEquals(b3, i3);
		input.close();
		f.delete();
	}

	public void testWriteUInt8() throws FileNotFoundException, IOException, UnsignedOutOfLimitsException {
		File f = new File("./TestFile.txt");
		BufferedFileOutputStream output = new BufferedFileOutputStream(f);
		UByte i1 = new UByte((short)67);
		FileOutput.writeUInt8(output, i1);
		UByte i2 = new UByte((short)96);
		FileOutput.writeUInt8(output, i2);
		UByte i3 = new UByte((short)54);
		FileOutput.writeUInt8(output, i3);
		output.close();

		FileInputStream input = new FileInputStream(f);
		UByte b1 = FileInput.readUInt8(input);
		assertEquals(b1.toShort(), i1.toShort());
		UByte b2 = FileInput.readUInt8(input);
		assertEquals(b2.toShort(), i2.toShort());
		UByte b3 = FileInput.readUInt8(input);
		assertEquals(b3.toShort(), i3.toShort());
		input.close();
		f.delete();
	}

	public void testWriteFloat() throws FileNotFoundException, IOException, UnsignedOutOfLimitsException {
		File f = new File("./TestFile.txt");
		BufferedFileOutputStream output = new BufferedFileOutputStream(f);
		float i1 = 2334.567f;
		FileOutput.writeFloat(output, i1);
		float i2 = 5778.96f;
		FileOutput.writeFloat(output, i2);
		float i3 = 37896.54f;
		FileOutput.writeFloat(output, i3);
		output.close();

		FileInputStream input = new FileInputStream(f);
		float b1 = FileInput.readFloat(input);
		assertEquals(b1, i1);
		float b2 = FileInput.readFloat(input);
		assertEquals(b2, i2);
		float b3 = FileInput.readFloat(input);
		assertEquals(b3, i3);
		input.close();
		f.delete();
	}

}
