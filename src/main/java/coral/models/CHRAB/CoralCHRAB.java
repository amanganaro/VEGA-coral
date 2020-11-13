/**
 * QSAR model for the chromosome aberration test.
 */
package coral.models.CHRAB;

import coral.CoralAttributes;
import coral.CoralModel;
import coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 * @author GLavado
 */
public class CoralCHRAB extends CoralModel{
    
    //Endpoint =  -0.3129603 (± 0.0019980) +    0.0682837 (± 0.0001541) * DCW(1,13)
    public CoralCHRAB() {
        super(-0.3129603, 0.0682837, "/attr-cw-CHRAB.txt");
    }

    @Override
    protected double calculateDCW(String smi) throws Exception {
        double DCW = 0;
        
        // list of tokens from the SMILES in input
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> list_tokens = coral_smi.ProcessSMILES(smi);
        
        // list of attributes
        ArrayList<String> list_attributes = new ArrayList<>();
       
        // Sk, SSk, HARD

        // build attributes Sk from the list of tokens
        list_attributes.addAll(CoralAttributes.calculateAttributeSks(list_tokens));
       
        // build attributes SSk from the list of tokens
        list_attributes.addAll(CoralAttributes.calculateAttributeSSks(list_tokens));
       
        list_attributes.add(CoralAttributes.calculateAttributeHARDCoral2017(list_tokens));
       
        // retrieve the weights of all attributes for the current loaded model
        ArrayList<Double> weights = calculateSAvalues(list_attributes);
       
       // ADDITIVE scheme for the optimal descriptor
        for (Double d : weights)
            DCW += d;
        
         return DCW;
    }
    
}
