package affymetrix.calvin.data;

import affymetrix.portability.DataSizes;

public class UIntColumn extends ColumnInfo {
	public UIntColumn(String name) {
		super(name, DataSetColumnTypes.UIntColType, DataSizes.INT_SIZE, 1, 0);
	}
}
