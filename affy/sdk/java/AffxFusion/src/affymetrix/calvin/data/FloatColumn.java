package affymetrix.calvin.data;

import affymetrix.portability.DataSizes;

public class FloatColumn extends ColumnInfo {
	public FloatColumn(String name) {
		super(name, DataSetColumnTypes.FloatColType, DataSizes.FLOAT_SIZE, 1, 0);
	}
}
