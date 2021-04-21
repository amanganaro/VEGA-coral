/*
* A Coral model for 204 compounds with pEC3 data was implemented under the project VERMEER.
* The experimental data is expressed in concentration (weight percent) giving a stimulation index value of 3 in the local lymph node assay (pEC3)  
* (as negative decimal logarithmic units). 
* The model involves Sk, SSk, SSSk and HARD attributes.
*/
package insilico.coral.models.skin.sensitization;


import insilico.coral.CoralAttributes;
import insilico.coral.CoralModel;
import insilico.coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 * @author Giovanna
 */
public class CoralSkinSensitization  extends CoralModel {

    public CoralSkinSensitization() {
        super(-0.4448002, 0.0860852, "/str_cw_sens.txt");
    }
    
    @Override
    protected double calculateDCW(String smi) throws Exception {
        // calculate the list of all Sk, SSk, SSSk and HARD attributes for the given SMILES
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> tokenList =  coral_smi.ProcessSMILES(smi);

        // Sk, SSk, SSSk
        ArrayList<String> attrList = CoralAttributes.calculateAllLocalAttributes(tokenList);
        
        attrList.add(CoralAttributes.calculateAttributeHARDCoral2017(tokenList));

        // retrieve the weights of all attributes for the current loaded model
        ArrayList<Double> weights = calculateSAvalues(attrList);
        
        // calculate the prediction with an ADDITIVE scheme
        double value = 0;
        for (Double d : weights)
            value += d;
        
        return value;
    }
    
}
