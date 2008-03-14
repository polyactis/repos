package affymetrix.calvin.data;

import affymetrix.portability.DataSizes;

public class UShortColumn extends ColumnInfo {
	public UShortColumn(String name) {
		super(name, DataSetColumnTypes.UShortColType, DataSizes.SHORT_SIZE, 1, 0);
	}
}
