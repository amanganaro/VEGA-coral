package coral.models.hydrolysis;

import coral.CoralAttributes;
import coral.CoralGraphAttributes;
import coral.CoralModel;
import coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 * @author GLavado
 */
public class CoralHydrolysis extends CoralModel{

    public CoralHydrolysis() {
        //Endpoint =  -1.2866535 (± 0.0418619) +    0.2544880 (± 0.0033488) * DCW(1,35)
        super(-1.2866535, 0.2544880, "/attr_cw_hydrolysis.txt");
    }

    @Override
    protected double calculateDCW(String smi) throws Exception {
       double DCW = 0;
       
        // Hydrogen suppressed graph (HSG) is used in this model
        // Attributes involve from graph: C6 C5
        // Attributes involve from SMILES: S(k), SS(k), SSS(k), (NO), (NS), (OS)
       
       // list of tokens from the SMILES in input
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> list_tokens = coral_smi.ProcessSMILES(smi);
        
        // sk, ssk, sssk
        ArrayList<String> list_attributes = CoralAttributes.calculateAllLocalAttributes(list_tokens);
        
        // (NO), (NS), (OS)
        String attr = CoralAttributes.calculateAttributesPairProportions(list_tokens, "N", "O");
        list_attributes.add(attr);
        attr = CoralAttributes.calculateAttributesPairProportions(list_tokens, "N", "S");
        list_attributes.add(attr);
        attr = CoralAttributes.calculateAttributesPairProportions(list_tokens, "O", "S");
        list_attributes.add(attr);
                
        // C6, C5
        CoralGraphAttributes coral_graph_attr = new CoralGraphAttributes(smi);
        
        for (int i=6; i>4;i--){
                String AttributeRings = coral_graph_attr.CalculateRings(i);
                list_attributes.add(AttributeRings);
            }
        /*
        for (String str: list_attributes){
            System.out.println(str);
        }
        */
         ArrayList<Double> weights = calculateSAvalues(list_attributes);
        
        // ADDITIVE scheme for the optimal descriptor
        for (Double d : weights){
            DCW += d;
        }
        //System.out.println(DCW);
        return DCW;
    }
    
}
