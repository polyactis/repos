import affymetrix.fusion.chp.FusionCHPData;
import affymetrix.fusion.chp.FusionCHPDataReg;
import affymetrix.fusion.chp.FusionCHPHeader;
import affymetrix.fusion.chp.FusionCHPLegacyData;
import affymetrix.fusion.chp.FusionExpressionProbeSetResults;

public class read_chp {

	public read_chp() {
	}

	/**
	 * @param args
	 *          the command line arguments
	 */
	public static void main(String[] args) {

		FusionCHPLegacyData.registerReader();
		String fileName = args[0];
		try {
			FusionCHPData chp = FusionCHPDataReg.read(fileName);
			if (chp == null) {
				System.out.println("Failed to read the file.");
				return;
			}
			FusionCHPLegacyData legchp = FusionCHPLegacyData.fromBase(chp);
			if (legchp == null) {
				System.out.println("The example is for expression CHP files only.");
				return;
			}
			int atype = legchp.getHeader().getAssayType();
			if (atype != FusionCHPHeader.EXPRESSION_ASSAY) {
				System.out.println("The example is for expression only.");
				return;
			}

			float sum = 0;
			int nump = 0;
			int n = legchp.getHeader().getNumProbeSets();
			FusionExpressionProbeSetResults psResults = new FusionExpressionProbeSetResults();
			for (int i = 0; i < n; i++) {
				legchp.getExpressionResults(i, psResults);
				sum += psResults.getSignal();
				if (psResults.getDetection().toShort() == FusionExpressionProbeSetResults.ABS_PRESENT_CALL) {
					++nump;
				}
			}
			float avg = sum / n;
			System.out.println("The average signal is " + avg);
			System.out.println("The #P is " + nump + " (" + (int)(100.0 * nump / n) + "%)");
		}
		catch (Exception e) {
		}
	}

}
