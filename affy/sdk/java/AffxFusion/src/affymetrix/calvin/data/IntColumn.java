package affymetrix.calvin.data;

import affymetrix.portability.DataSizes;

public class IntColumn extends ColumnInfo {
	public IntColumn(String name_) {
		super(name_, DataSetColumnTypes.IntColType, DataSizes.INT_SIZE, 1, 0);
	}
}
