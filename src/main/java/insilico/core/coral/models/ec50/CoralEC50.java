/**
 * CORAL model for 252 compounds with experimental data for general NOAEL in Rats for oral sub-chronic RDT. 
 * The experimental data is expressed in mg/kg bw (as logarithmic units)
 * The optimal descriptor is based on attributes extracted from SMILES (Sk, SSk, SSSk, e1, e2, e3, nn).
 */
package insilico.core.coral.models.ec50;

import insilico.core.coral.CoralAttributes;
import insilico.core.coral.CoralModel;
import insilico.core.coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 */
public class CoralEC50 extends CoralModel{


    public CoralEC50() {
        super(-0.8972636 ,  0.1447300, "/attr-cw-ec50-corr.txt");
    }
    
        @Override
    protected double calculateDCW(String smi) throws Exception {
        double dcw = 0;

        // tokens from SMILES string
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> res = coral_smi.ProcessSMILESwithSquareBracket(smi);

        // Sk,SSk,SSSk on all elements NOT including square brackets
        CoralAttributes.SQUARED_BRACKETS_TOKENS_CALCULATION = false;
        ArrayList<String> attr_list = CoralAttributes.calculateAllLocalAttributes(res);

        ArrayList<Double> weights = calculateSAvalues(attr_list, false);
        for (Double cw : weights) {
            dcw += cw;
        }
        return dcw;

    }
}
    

