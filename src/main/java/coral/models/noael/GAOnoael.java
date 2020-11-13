/**
 * A Coral model for 140 compounds with experimental data on NOAEL (expressed as -Log mg/kg)
 * The optimal descriptor is based on attributes extracted from SMILES (BOND, NOSP, HALO) and GAO descriptors (EC0).
 * Paper: CORAL: model for no observed adverse effect level (NOAEL)”, Toropov et al, 2015.
 * 
 */
package coral.models.noael;

import coral.CoralAttributes;
import coral.CoralModel;
import coral.CoralSMILES;
import coral.GAODescriptor;

import java.util.ArrayList;

/**
 *
 * @author GLavado
 */
public class GAOnoael extends CoralModel {

    public GAOnoael() {
        super(-2.1680835, 0.0737528, "/GAO_noael_attr_cw.txt");
    }

    @Override
    protected double calculateDCW(String smi) throws Exception {
        
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> tokenList =  coral_smi.ProcessSMILES(smi);
        
        GAODescriptor gao = new GAODescriptor();
        gao.SetHAddType(GAODescriptor.H_TYPE_AROMATIC);
        
        // calculate adjacency matrix GAO
        gao.Calculate(smi);
        
        // calculate EC0 from the given GAO
        ArrayList<String> attrList = gao.CalculateEC0();
           
        // calculate NOSP, HALO, BOND attributes for the given SMILES
        attrList.add(CoralAttributes.calculateAttributeNOSPCoral2017(tokenList));
        attrList.add(CoralAttributes.calculateAttributeHALOCoral2017(tokenList));
        attrList.add(CoralAttributes.calculateAttributeBOND(tokenList));
        
        /*for (String att : attrList) {
            System.out.println(att);
        }
        */
        
        // retrieve the weights of all attributes for the current loaded model
        ArrayList<Double> weights = calculateSAvalues(attrList);
        
        // calculate the prediction with an ADDITIVE scheme
        double value = 0;
        for (Double d : weights)
            value += d;

        return value;
    }
    
}
