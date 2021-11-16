
package insilico.core.coral.models.airhalflife;

import insilico.core.coral.CoralAttributes;
import insilico.core.coral.CoralModel;
import insilico.core.coral.CoralSMILES;

import java.util.ArrayList;

/**
 * A Coral model for 302 compounds with experimental data on half-life for Air (expressed as Log hours in base 10)
 * was developed under the project VERMEER. The optimal descriptor is based on attributes extracted from SMILES (Sk, SSk, SSSk, HARD).
 *
 */
public class CoralAirHalfLife extends CoralModel {

    public CoralAirHalfLife() {
        super(-2.2271129,  0.0900047, "/list_airhalflife_SA_CW.txt");
    }

    @Override
    protected double calculateDCW(String smi) throws Exception {
        double DCW = 0; 
        
        // list of tokens from the SMILES in input
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> list_tokens = coral_smi.ProcessSMILES(smi);
        
        // Sk, SSk and SSSk attributes from SMILES
        ArrayList<String> list_attributes = CoralAttributes.calculateAttributeSks(list_tokens);
        list_attributes.addAll(CoralAttributes.calculateAttributeSSks(list_tokens));
        list_attributes.addAll(CoralAttributes.calculateAttributeSSSks(list_tokens));
        
        list_attributes.add(CoralAttributes.calculateAttributeHARDCoral2017(list_tokens));
        
        // retrieve the weights of all attributes for the current loaded model
        ArrayList<Double> weights = calculateSAvalues(list_attributes);
        
        // ADDITIVE scheme for the optimal descriptor
        for (Double d : weights){
            DCW += d;
        }
            
        return DCW;
    }

    
}
