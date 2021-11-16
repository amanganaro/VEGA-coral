package insilico.core.coral.models.noael;

import insilico.core.coral.CoralAttributes;
import insilico.core.coral.CoralModel;
import insilico.core.coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 * @author GLavado
 */
public class CoralNOAEL extends CoralModel {

    public CoralNOAEL() {
        super(1.1574991, 0.0608755 , "/noael_attr_cw.txt");
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
