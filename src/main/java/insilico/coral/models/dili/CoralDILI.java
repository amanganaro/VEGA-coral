/*
 * A Coral model for human hepatotoxicity (drug-induced liver injury) was implemented under the project VERMEER. 
 * The dataset consists of 2029 therapeutic agents with the drug-induced liver injury (DILI) expressed 
 * as active (1) and inactive (0).

 * The model involves HSG and SMILES attributes: EC0, EC1, EC2, NNC, C3, C5, C6, Sk, SSk, SSSk, HARD.
 *
 * Paper: 
 * Alla P. Toropova, Andrey A. Toropov,
 * CORAL: Binary classifications (active/inactive) for drug-induced liver injury, Toxicology Letters, Volume 268, 2017, Pages 51-57, ISSN 0378-4274
 *
 */

package insilico.coral.models.dili;


import insilico.coral.CoralAttributes;
import insilico.coral.CoralGraphAttributes;
import insilico.coral.CoralModel;
import insilico.coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 * @author Giovanna
 */
public class CoralDILI extends CoralModel {

    public CoralDILI() {
        // Endpoint =  -0.1390926 (± 0.0006896) +    0.0163083 (± 0.0000188) * DCW(1,30)
        super(-0.1390926,  0.0163083, "/dili/liver_str-cw.txt");
    }

    @Override
    protected double calculateDCW(String smi) throws Exception {

        // HSG attributes: EC0, EC1, EC2, NNC, C3, C5, C6
       CoralGraphAttributes coral_graph_attr = new CoralGraphAttributes(smi);
        
       // SMILES attributes: Sk, SSk, SSSk, HARD
       CoralSMILES coral_smi = new CoralSMILES();
       ArrayList<String> tokenList =  coral_smi.ProcessSMILES(smi);
        
       ArrayList<String> attrList = new ArrayList<>();
 
        // EC0, EC1, EC2
       attrList.addAll(coral_graph_attr.CalculateEC(2));
        // NNC
       attrList.addAll(coral_graph_attr.CalculateNNC());
       
       String attributeRings;
       // C6, C5
       for (int i=6; i>4;i--){
           //attributeRings = coral_graph_attr.CalculateRings(i);
           attributeRings = coral_graph_attr.CalculateRingsWithHeteroAtoms(i);
           attrList.add(attributeRings);
       }
       // C3
       //attributeRings = coral_graph_attr.CalculateRings(3);
       attributeRings = coral_graph_attr.CalculateRingsWithHeteroAtoms(3);
       attrList.add(attributeRings);
        
       // calculate the list of all Sk, SSk, SSSk attributes for the given SMILES
       attrList.addAll(CoralAttributes.calculateAllLocalAttributes(tokenList));
        
       // calculate HARD attribute from SMILES
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
