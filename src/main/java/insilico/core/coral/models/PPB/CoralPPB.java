/**
 * A Coral model for 489 compounds with experimental data for plasma protein binding (PPB) was developed under the project EU-ToxRisk. 
 * The experimental data is expressed as square root of fraction unbound (fu)
 * The optimal descriptor is based on attributes extracted from SMILES (Sk, SSk, SSSk, HARD).
 * Paper:
 * Toma, C., Gadaleta, D., Roncaglioni, A. et al. 
 * QSAR Development for Plasma Protein Binding: Influence of the Ionization State. Pharm Res 36, 28 (2019). 
 */
package insilico.core.coral.models.PPB;

import insilico.core.coral.CoralAttributes;
import insilico.core.coral.CoralModel;
import insilico.core.coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 * @author GLavado
 */
public class CoralPPB extends CoralModel {

    // Endpoint =   0.7454093 (± 0.0017160) +    0.0213199 (± 0.0000984) * DCW(1,1)
    public CoralPPB() {
        super(0.7454093, 0.0213199, "/ppb_str_cw.txt");
    }

    @Override
    protected double calculateDCW(String smi) throws Exception {
       
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> tokenList =  coral_smi.ProcessSMILES(smi);
        
        // calculate the list of all Sk, SSk and SSSk for the given SMILES
        ArrayList<String> attrList = CoralAttributes.calculateAllLocalAttributes(tokenList);
        // HARD attribute
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
