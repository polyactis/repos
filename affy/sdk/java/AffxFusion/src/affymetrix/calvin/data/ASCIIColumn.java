package affymetrix.calvin.data;

import affymetrix.portability.DataSizes;

public class ASCIIColumn extends ColumnInfo {
	public ASCIIColumn(String name, int maxLn) {
		super(name, DataSetColumnTypes.ASCIICharColType, DataSizes.CHAR_SIZE, maxLn, 4);
	}
}
