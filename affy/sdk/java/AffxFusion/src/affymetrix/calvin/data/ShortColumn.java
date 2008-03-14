package affymetrix.calvin.data;

import affymetrix.portability.DataSizes;

public class ShortColumn extends ColumnInfo {
	public ShortColumn(String name) {
		super(name, DataSetColumnTypes.ShortColType, DataSizes.SHORT_SIZE, 1, 0);
	}
}
