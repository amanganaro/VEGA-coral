/**
 * CORAL model for 252 compounds with experimental data for general NOAEL in Rats for oral sub-chronic RDT. 
 * The experimental data is expressed in mg/kg bw (as logarithmic units)
 * The optimal descriptor is based on attributes extracted from SMILES (Sk, SSk, SSSk, e1, e2, e3, nn).
 */
package insilico.core.coral.models.noael.liver;

import insilico.core.coral.CoralAttributes;
import insilico.core.coral.CoralGraphAttributes;
import insilico.core.coral.CoralModel;
import insilico.core.coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 * @author MMarzo
 */
public class CoralLiverNoael extends CoralModel{

    //Endpoint =   1.7245212 (± 0.0040067) +    0.0530130 (± 0.0002568) * DCW(2,15)
    public CoralLiverNoael() {
        super(1.7245212, 0.0530130, "/attr-cw-noael-liver.txt");
    }
    
        @Override
    protected double calculateDCW(String smi) throws Exception {
        double dcw = 0;
        
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> list_tokens = coral_smi.ProcessSMILES(smi);
        
        
        //Sk, SSk, SSSk, e1, e2, e3, nn
        ArrayList<String> list_attributes = CoralAttributes.calculateAllLocalAttributes(list_tokens);
        
      
        CoralGraphAttributes coral_graph_attr = new CoralGraphAttributes(smi);
        // NNC
        list_attributes.addAll(coral_graph_attr.CalculateNNC());
        // EC1, EC2, EC3
        list_attributes.addAll(coral_graph_attr.CalculateEC(3));
        
        ArrayList<Double> weights = calculateSAvalues(list_attributes);
        
        // ADDITIVE scheme for the optimal descriptor
        for (Double d : weights){
            dcw += d;
        }
        return dcw;
    }
}
    

