package be.iminds.cloudspeller.output;

import be.iminds.cloudspeller.motifmodels.FreqVec;
import be.iminds.cloudspeller.phylogenetics.BLS;
import be.iminds.cloudspeller.toolbox.GeneralToolbox;

/**
 * Created by ddewitte on 09.07.15.
 */
public class TripleThresholdFilterExactBLST  extends SimultaneousOccurrenceConfidenceAndBLSFiltering {


    public TripleThresholdFilterExactBLST(int confidenceCutoff, int familyOccurrenceCutoff, int blsCutoff) {
        super(confidenceCutoff, familyOccurrenceCutoff, blsCutoff);
    }


    @Override
    public boolean checkRestrictions(String strValue) {

        GeneralToolbox.parseConfidenceGraphValues(strValue, F, p);
        int i=minBLSId;
        if (p[i]>=probCutoff){
                if (F[i]>=famCutoff)
                    return true;

        }
        return false;

    }

    @Override
    public boolean checkRestrictions(int [] F, double [] p){
        this.F = F;
        this.p = p;
        int i=minBLSId;
            if (p[i]>=probCutoff){
                if (F[i]>=famCutoff)
                    return true;

        }
        return false;

    }

}
