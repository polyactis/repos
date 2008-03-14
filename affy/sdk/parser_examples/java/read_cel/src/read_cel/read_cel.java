import affymetrix.fusion.cdf.FusionCDFData;
import affymetrix.fusion.cdf.FusionCDFProbeGroupInformation;
import affymetrix.fusion.cdf.FusionCDFProbeInformation;
import affymetrix.fusion.cdf.FusionCDFProbeSetInformation;
import affymetrix.fusion.cel.FusionCELData;

public class read_cel {

	public read_cel() {
	}

	/**
	 * @param args
	 *          the command line arguments
	 */
	public static void main(String[] args) {

		String celFileName = args[0];
		String cdfFileName = args[1];
		FusionCELData cel;
		FusionCDFData cdf;
		try {
			cel = new FusionCELData();
			cel.setFileName(celFileName);
			if (cel.read() == false) {
				System.out.println("Failed to read the CEL file.");
				return;
			}

			float sum = 0;
			int n = cel.getCells();
			for (int i = 0; i < n; i++) {
				sum += cel.getIntensity(i);
			}
			float avg = sum / n;
			System.out.println("The average intensity is " + avg);

			cdf = new FusionCDFData();
			cdf.setFileName(cdfFileName);
			if (cdf.read() == false) {
				System.out.println("Failed to read the CDF file.");
				return;
			}

			int nsets = cdf.getHeader().getNumProbeSets();
			for (int iset = 0; iset < nsets; iset++) {
				sum = 0;
				FusionCDFProbeSetInformation set = new FusionCDFProbeSetInformation();
				cdf.getProbeSetInformation(iset, set);
				int ngroups = set.getNumGroups();
				for (int igroup = 0; igroup < ngroups; igroup++) {
					FusionCDFProbeGroupInformation group = new FusionCDFProbeGroupInformation();
					set.getGroup(igroup, group);
					int ncells = group.getNumCells();
					for (int icell = 0; icell < ncells; icell++) {
						FusionCDFProbeInformation probe = new FusionCDFProbeInformation();
						group.getCell(icell, probe);
						sum += cel.getIntensity(probe.getX(), probe.getY());
					}
				}
				avg = sum / set.getNumCells();
				System.out.println("The average probe set intensity is " + avg);
			}
		}
		catch (Exception e) {
		}
	}
}
