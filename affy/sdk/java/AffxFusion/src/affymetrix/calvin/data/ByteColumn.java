package affymetrix.calvin.data;

import affymetrix.portability.DataSizes;

public class ByteColumn extends ColumnInfo {
	public ByteColumn(String name) {
		super(name, DataSetColumnTypes.ByteColType, DataSizes.CHAR_SIZE, 1, 0);
	}
}
