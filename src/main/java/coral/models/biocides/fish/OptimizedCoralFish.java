/**
 * A Coral model for 88 biocides with experimental data on LC50 for Fish Acute Toxicity (expressed as –Log mmol/L in base 10)
 * was developed under the project COMBASE. The optimal descriptor is based on local attributes (Sk, SSk, SSSk).
 * 
 */

package insilico.core.coral.models.biocides.fish;

import insilico.core.coral.CoralModel;
import insilico.core.coral.CoralSMILES;
import java.util.ArrayList;
import insilico.core.coral.CoralAttributes;
/**
 *
 * @author Giovanna
 */
public class OptimizedCoralFish extends CoralModel{

    public OptimizedCoralFish() {
        super(-0.9887356,  0.0454373, "/insilico/core/coral/models/biocides/fish/optimized_fish_attr_cw.txt");
    }

    @Override
    protected double calculateDCW(String smi) throws Exception {
        
        CoralSMILES coral_smiles = new CoralSMILES();
        ArrayList<String> list_tokens =  coral_smiles.ProcessSMILES(smi);
        
        // Sk, SSk and SSSk attributes from SMILES
        ArrayList<String> list_attributes = CoralAttributes.calculateAttributeSks(list_tokens);
        list_attributes.addAll(CoralAttributes.calculateAttributeSSks(list_tokens));
        list_attributes.addAll(CoralAttributes.calculateAttributeSSSks(list_tokens));
            
        // retrieve the weights of all attributes for the current loaded model
        ArrayList<Double> weights = calculateSAvalues(list_attributes);
        
        // ADDITIVE scheme for the optimal descriptor
        double dcw = 0;
        for (Double d : weights)
            dcw += d;
        
        return dcw;
    }
    
}
