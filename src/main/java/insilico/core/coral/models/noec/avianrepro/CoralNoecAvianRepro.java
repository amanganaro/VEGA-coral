/**
 * CORAL model for 252 compounds with experimental data for general NOAEL in Rats for oral sub-chronic RDT. 
 * The experimental data is expressed in mg/kg bw (as logarithmic units)
 * The optimal descriptor is based on attributes extracted from SMILES (Sk, SSk, SSSk, e1, e2, e3, nn).
 */
package insilico.core.coral.models.noec.avianrepro;

import insilico.core.coral.CoralAttributes;
import insilico.core.coral.CoralModel;
import insilico.core.coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 */
public class CoralNoecAvianRepro extends CoralModel{


    public CoralNoecAvianRepro() {
        super(-1.2755261 ,  0.0572993, "/attr-cw-noec-avian-repro.txt");
    }
    
        @Override
    protected double calculateDCW(String smi) throws Exception {
        double dcw = 0;

        // tokens from SMILES string
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> res = coral_smi.ProcessSMILESwithSquareBracket(smi);

        // Sk,SSk,SSSk on all elements NOT including square brackets
        CoralAttributes.SQUARED_BRACKETS_TOKENS_CALCULATION = true;
//        ArrayList<String> attr_list = CoralAttributes.calculateAllLocalAttributes(res);
        ArrayList<String> attr_list = CoralAttributes.calculateAttributeSks(res);
        attr_list.addAll( CoralAttributes.calculateAttributeSSks(res) );
        // symmetry attributes are computed on the raw SMILES characters, NOT on the tokens
        // (see CoralAttributes.calculateAttributeSymmetry for the details)
        attr_list.addAll(CoralAttributes.calculateAttributeSymmetry(smi));
        CoralAttributes.SQUARED_BRACKETS_TOKENS_CALCULATION = false;

        ArrayList<Double> weights = calculateSAvalues(attr_list, false);
        for (Double cw : weights) {
            dcw += cw;
        }
        return dcw;

    }
}
    

