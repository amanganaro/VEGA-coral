/**
 * CORAL model for 252 compounds with experimental data for general NOAEL in Rats for oral sub-chronic RDT. 
 * The experimental data is expressed in mg/kg bw (as logarithmic units)
 * The optimal descriptor is based on attributes extracted from SMILES (Sk, SSk, SSSk, e1, e2, e3, nn).
 */
package insilico.core.coral.models.henryslaw;

import insilico.core.coral.CoralAttributes;
import insilico.core.coral.CoralGraphAttributes;
import insilico.core.coral.CoralModel;
import insilico.core.coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 * @author MMarzo
 */
public class CoralHenrysLaw extends CoralModel{


    public CoralHenrysLaw() {
        super(3.8465029,  1.1055325, "/attr-cw-henry_law.txt");
    }
    
        @Override
    protected double calculateDCW(String smi) throws Exception {
        double dcw = 0;

        // tokens from SMILES string
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> res = coral_smi.ProcessSMILESwithSquareBracket(smi);

        // Sk,SSk,SSSk on all elements NOT including square brackets
        CoralAttributes.SQUARED_BRACKETS_TOKENS_CALCULATION = true;
        ArrayList<String> attr_list = CoralAttributes.calculateAllLocalAttributes(res);
        CoralAttributes.SQUARED_BRACKETS_TOKENS_CALCULATION = false;



        ArrayList<Double> weights = calculateSAvalues(attr_list, false);
        for (Double cw : weights) {
            dcw += cw;
        }
        return dcw;

    }
}
    

