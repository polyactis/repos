

package fusionsample;

import affymetrix.fusion.psi.*;
import affymetrix.fusion.cel.*;
import affymetrix.fusion.chp.*;
import affymetrix.fusion.*;
import java.util.*;

/** This class provides sample code for using the Fusion SDK.
 * Examples provided include:
 * 1. exportPSI - Reads the PSI file and displays the probe set names.
 * 2. exportCHP - Reads the CHP file for expression and genotyping only and outputs some summary statistics.
 * 3. exportCEL - Reads the CEL file and outputs the average intensity value.
 */
public class Main {
    
    /** Creates a new instance of Main */
    public Main() {
    }
    
    /** Calculates some statistics (# and % calls and average signal value) for expression probe set results.
     * @param chp The CHP file object.
     */
    private void exportExpressionCHP(FusionCHPLegacyData chp) {
        int numPresent=0;
        int numAbsent=0;
        int numMarginal=0;
        double signalPresent=0.0;
        double signalAbsent=0.0;
        double signalMarginal=0.0;
        int n = chp.getHeader().getNumProbeSets();
        FusionExpressionProbeSetResults result = new FusionExpressionProbeSetResults();
        for (int i=0; i<n; i++)
        {
            chp.getExpressionResults(i, result);
            switch (result.getDetection())
            {
                case FusionExpressionProbeSetResults.ABS_ABSENT_CALL:
                    ++numAbsent;
                    signalAbsent += result.getSignal();
                    break;

                case FusionExpressionProbeSetResults.ABS_PRESENT_CALL:
                    ++numPresent;
                    signalPresent += result.getSignal();
                    break;

                case FusionExpressionProbeSetResults.ABS_MARGINAL_CALL:
                    ++numMarginal;
                    signalMarginal += result.getSignal();
                    break;

                default:
                    break;
            }
        }
        System.out.println("Present calls:\t" + numPresent + "\t" + (100.0f*numPresent/n) + "%");
        System.out.println("Present signal:\t" + signalPresent/numPresent);
        System.out.println("Absent calls:\t" + numAbsent + "\t" + (100.0f*numAbsent/n) + "%");
        System.out.println("Absent signal:\t" + signalAbsent/numAbsent);
        System.out.println("Marginal calls:\t" + numMarginal + "\t" + (100.0f*numMarginal/n) + "%");
        System.out.println("Marginal signal:\t" + signalMarginal/numMarginal);
    }

    /** Calculates some statistics (# and % calls) for genotyping probe set results.
     * @param chp The CHP file object.
     */
    private void exportGenotypingCHP(FusionCHPLegacyData chp) {
        int numAA=0;
        int numAB=0;
        int numBB=0;
        int n = chp.getHeader().getNumProbeSets();
        FusionGenotypeProbeSetResults result = new FusionGenotypeProbeSetResults();
        for (int i=0; i<n; i++)
        {
            chp.getGenotypingResults(i, result);
            switch (result.getAlleleCall())
            {
                case FusionGenotypeProbeSetResults.ALLELE_A_CALL:
                    ++numAA;
                    break;
                    
                case FusionGenotypeProbeSetResults.ALLELE_B_CALL:
                    ++numBB;
                    break;
                case FusionGenotypeProbeSetResults.ALLELE_AB_CALL:
                    ++numAB;
                    break;
                    
                default:
            }
        }
        System.out.println("AA calls:\t" + numAA + "\t" + (100.0f*numAA/n) + "%");
        System.out.println("AB calls:\t" + numAB + "\t" + (100.0f*numAB/n) + "%");
        System.out.println("BB calls:\t" + numBB + "\t" + (100.0f*numBB/n) + "%");
    }
    
    /** Exports data stored in a CHP file.
     * @param file The name of the CHP file.
     */
    private void exportCHP(String file) {
        FusionCHPLegacyData.registerReader();
        FusionCHPData chp = FusionCHPDataReg.read(file);
        if (chp == null)
        {
            System.out.println("Failed to read the CHP file.");
            return;
        }        
        FusionCHPLegacyData fusionchp = FusionCHPLegacyData.fromBase(chp);
        
        System.out.println("Alg name = " + fusionchp.getHeader().getAlgName());
        System.out.println("Alg version = " + fusionchp.getHeader().getAlgVersion());
        System.out.println("Array type = " + fusionchp.getHeader().getChipType());
        System.out.println("# probe sets = " + fusionchp.getHeader().getNumProbeSets());

        Vector params = fusionchp.getHeader().getAlgorithmParameters();
        for (int i=0; i<params.size(); i++)
        {
            FusionTagValuePair param = (FusionTagValuePair) params.elementAt(i);
            System.out.println(param.getTag() + " = " + param.getValue());
        }
        params = fusionchp.getHeader().getSummaryParameters();
        for (int i=0; i<params.size(); i++)
        {
            FusionTagValuePair param = (FusionTagValuePair) params.elementAt(i);
            System.out.println(param.getTag() + " = " + param.getValue());
        }
                
        if (fusionchp.getHeader().getAssayType() == FusionCHPHeader.EXPRESSION_ASSAY)
        {
            exportExpressionCHP(fusionchp);
        }
        else if (fusionchp.getHeader().getAssayType() == FusionCHPHeader.GENOTYPING_ASSAY)
        {
            exportGenotypingCHP(fusionchp);
        }
    }
    
    /** Exports the probe set names stored in a PSI file.
     * @param file The PSI file name.
     * @param startIndex The 0 based index of the first probe set to print.
     * @param endIndex The 0 based index of the last probe set to print. -1 for all probe sets after the start index.
     */
    private void exportPSI(String file, int startIndex, int endIndex) {
        FusionPSIData psi = new FusionPSIData();
        psi.setFileName(file);
        if (psi.exists() == false)
        {
            System.out.println("The PSI file does not exist.");
            return;
        }
        if (psi.read() == false)
        {
            System.out.println("Failed to read the PSI file.");
            return;
        }
        int n = psi.getProbeSetCount();
        int start = Math.min(startIndex, n-1);
        int end = Math.min(endIndex+1,n);
        if (end == 0)
            end = n;
        System.out.println("There are " + n + " probe sets in the PSI file.");
        System.out.println("Showing probe set names " + (start+1) + " to " + end);
        for (int i=start; i<end; i++)
            System.out.println("Probe set #" + (i+1) + "\t" + psi.getProbeSetName(i));
    }
    
    /** Exports the average intensity in a CEL file.
     * @param file The name of the CEL file.
     */
    private void exportCEL(String file) {
        FusionCELData cel = new FusionCELData();
        
        cel.setFileName(file);
        if (cel.exists() == false)
        {
            System.out.println("The CEL file does not exist.");
            return;
        }
        if (cel.read() == false)
        {
            System.out.println("Failed to read the CEL file.");
            return;
        }
        
        System.out.println("Algorithm name = " + cel.getAlg());
        System.out.println("Array type = " + cel.getChipType());
        Vector params = cel.getParameters();
        for (int i=0; i<params.size(); i++)
        {
            FusionTagValuePair param = (FusionTagValuePair) params.elementAt(i);
            System.out.println(param.getTag() + " = " + param.getValue());
        }
        
        int n = cel.getCells();
        double sum = 0;
        for (int i=0; i<n; i++)
            sum += cel.getIntensity(i);
        System.out.println("The average intensity is " + (int)(sum/n));
    }
    
    private String getFile(String[] args) {
        String file = null;
        for (int i=0; i<args.length; i++)
        {
            // Get the file name.
            if (args[i].compareToIgnoreCase("-file") == 0)
            {
                file = args[i+1];
                break;
            }
        }
        return file;
    }
    
    /** Calls the various export functions for PSI, CHP, CEL and CDF data.
     * @param args the command line arguments
     */
    public void exportFiles(String[] args) {
        // Find the input file name.
        String file = getFile(args);
        if (file == null)
            return;

        String ext = "psi";
        int toffset = file.length()-ext.length();
        if (file.regionMatches(true, toffset, ext, 0, ext.length()) == true)
        {
            int startIndex = 0;
            int endIndex = -1;
            exportPSI(file, startIndex, endIndex);
        }

        ext = "chp";
        toffset = file.length()-ext.length();
        if (file.regionMatches(true, toffset, ext, 0, ext.length()) == true)
        {
            exportCHP(file);
        }

        ext = "cel";
        toffset = file.length()-ext.length();
        if (file.regionMatches(true, toffset, ext, 0, ext.length()) == true)
        {
            exportCEL(file);
        }
    }

    /** Calls the various export functions for PSI, CHP, CEL and CDF data.
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        Main sample = new Main();
        sample.exportFiles(args);
    }    
}
