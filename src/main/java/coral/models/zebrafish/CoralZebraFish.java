package coral.models.zebrafish;

import coral.CoralAttributes;
import coral.CoralModel;
import coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 * @author Giovanna Lavado
 * @author Alberto Manganaro
 */
public class CoralZebraFish extends CoralModel {
    
    public CoralZebraFish() throws Exception {
        super(1.1053641, 0.0272434, "/zebrafish-2.txt");
    }

    
    @Override
    protected double calculateDCW(String smi) throws Exception {
        
        // calculate the list of all Sk, SSk and SSSk for the given SMILES
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> tokenList =  coral_smi.ProcessSMILES(smi);
        ArrayList<String> localAttrList = CoralAttributes.calculateAllLocalAttributes(tokenList);
    
        // retrieve the weights of all attributes for the current loaded model
        ArrayList<Double> weights = calculateSAvalues(localAttrList);
        
        // calculate the prediction with an ADDITIVE scheme
        double value = 0;
        for (Double d : weights)
            value += d;

        return value;        
    }
    
}
