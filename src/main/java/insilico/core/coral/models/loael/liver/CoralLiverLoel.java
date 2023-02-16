/**
 * CORAL model for 252 compounds with experimental data for liver LOAEL in Rats for oral sub-chronic RDT. 
 * The experimental data is expressed in mg/kg bw (as logarithmic units)
 * The optimal descriptor is based on attributes extracted from SMILES (Sk, SSk, SSSk, e1, e2, e3, nn).
 */

package insilico.core.coral.models.loael.liver;

import insilico.core.coral.*;

import java.util.ArrayList;

/**
 *
 * @author MIRKO-PC
 */
public class CoralLiverLoel extends CoralModel{
    
    //Endpoint =   1.7148691 (± 0.0038133) +    0.0529304 (± 0.0002685) * DCW(2,14)
    public CoralLiverLoel() {
        super (1.7148691,0.0529304,"/attr-cw-loael-liver.txt");
    }

    @Override
    protected double calculateDCW(String smi) throws Exception {
        double dcw = 0;

        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> list_tokens = coral_smi.ProcessSMILES(smi);

        // Sk, SSk, SSSk, e1, e2, e3, nn
        ArrayList<String> list_attributes = CoralAttributes.calculateAllLocalAttributes(list_tokens);

        CoralGraphAttributes coral_graph_attr = new CoralGraphAttributes(smi);
        // NNC
        list_attributes.addAll(coral_graph_attr.CalculateNNC());
        // EC,EC2,EC3
        list_attributes.addAll(coral_graph_attr.CalculateEC(3));

        ArrayList<Double> weights = calculateSAvalues(list_attributes);

        // ADDITIVE scheme for the optimal descriptor
        for (Double d : weights){
            dcw += d;
        }

        return dcw;

    }

}
