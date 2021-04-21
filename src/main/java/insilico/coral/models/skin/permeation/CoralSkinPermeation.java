/*
* A Coral model for a dataset of 271 compounds with Kp (skin permeability coefficient) data was implemented under the project VERMEER.
* The experimental data is expressed cm/h with a logarithm transformation.  
* The model involves Sk, SSk and SSSk attributes.
 */
package insilico.coral.models.skin.permeation;

import insilico.coral.CoralAttributes;
import insilico.coral.CoralModel;
import insilico.coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 * @author GLavado
 */
public class CoralSkinPermeation extends CoralModel {

    public CoralSkinPermeation() {
        // intercept, slope, SaSource
        super(-2.7829828, 0.2139683, "/skin_perm_attr_cw.txt");
    }
    
    @Override
    protected double calculateDCW(String smi) throws Exception {
        // calculate the list of all Sk, SSk and SSSk attributes for the given SMILES
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> tokenList =  coral_smi.ProcessSMILES(smi);

        // Sk, SSk, SSSk
        ArrayList<String> attrList = CoralAttributes.calculateAllLocalAttributes(tokenList);
        
        // retrieve the weights of all attributes for the current loaded model
        ArrayList<Double> weights = calculateSAvalues(attrList);
        
        // calculate the prediction with an ADDITIVE scheme
        double value = 0;
        for (Double d : weights)
            value += d;
        
        return value;
    }
    
}
