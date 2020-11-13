/**
 * The Coral model for 133 biocides with experimental data on EC50 48h Daphnia magna (expressed as Log mmol/L)
 * was developed under the project COMBASE.  The model uses local attributes (Sk, SSk).
 * 
 */
package coral.models.biocides.daphniamagna;

import coral.CoralAttributes;
import coral.CoralModel;
import coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 * @author GLavado
 */
public class OptimizedCoralDaphnia extends CoralModel {

    public OptimizedCoralDaphnia() {
        // intercept, slope, structural attributes with their correlation weights
        super(0.5569253, 0.1098193 , "/optimized_daphnia_SA_CW.txt");
    }

    @Override
    protected double calculateDCW(String smi) throws Exception {
       double dcw = 0;
       
       // extract tokens from SMILES smi
       CoralSMILES coral_smiles = new CoralSMILES();
       ArrayList<String> list_tokens =  coral_smiles.ProcessSMILES(smi);
       
       // list of attributes
       ArrayList<String> list_attributes = new ArrayList<>();
       
       // build attributes Sk from the list of tokens
       list_attributes.addAll(CoralAttributes.calculateAttributeSks(list_tokens));
       
       // build attributes SSk from the list of tokens
       list_attributes.addAll(CoralAttributes.calculateAttributeSSks(list_tokens));
       
       // retrieve the weights of all attributes for the current loaded model
        ArrayList<Double> weights = calculateSAvalues(list_attributes);
       
       // ADDITIVE scheme for the optimal descriptor
        for (Double d : weights)
            dcw += d;
        
        return dcw;
    }
    
}
