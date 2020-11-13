package coral.models.loael;


import coral.CoralAttributes;
import coral.CoralModel;
import coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 * @author GLavado
 */
public class CoralLOAEL extends CoralModel{

    public CoralLOAEL() throws Exception {
        super(1.6712499, 0.0567754, "/loael_attr_cw.txt");
    }

    @Override
    protected double calculateDCW(String smi) throws Exception {
        
        // calculate the list of all Sk, SSk and HARD attributes for the given SMILES
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> tokenList =  coral_smi.ProcessSMILES(smi);
        
        ArrayList<String> attrList = CoralAttributes.calculateAttributeSks(tokenList);
        attrList.addAll(CoralAttributes.calculateAttributeSSks(tokenList));
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
