package affymetrix.calvin.data;

import affymetrix.portability.DataSizes;

public class UnicodeColumn extends ColumnInfo {
	public UnicodeColumn(String name, int maxLn) {
		super(name, DataSetColumnTypes.UnicodeCharColType, DataSizes.SHORT_SIZE, maxLn, 4);
	}
}
