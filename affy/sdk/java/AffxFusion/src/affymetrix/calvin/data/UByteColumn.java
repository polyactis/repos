package affymetrix.calvin.data;

import affymetrix.portability.DataSizes;

public class UByteColumn extends ColumnInfo {
	public UByteColumn(String name) {
		super(name, DataSetColumnTypes.UByteColType, DataSizes.CHAR_SIZE, 1, 0);
	}
}
