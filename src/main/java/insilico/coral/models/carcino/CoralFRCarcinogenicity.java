/*
 * A Coral model for carcinogenic potency expressed as tumorigenic dose, TD50, was implemented under the project CONCERT REACH.
 * The numerical data on TD50 are expressed as mg/kg bw/day with a negative decimal logarithm transformation, i.e. pTD50.
 * The female rat (FR) dataset consists of 158 compounds.
 * The model involves HSG and SMILES attributes: EC0, NNC, C5, C6, S, SS, BOND, NOSP,HALO.
 *
 * Paper: 
 * Toropova, A.P., Toropov, A.A. CORAL: QSAR models for carcinogenicity of organic compounds for male and female rats (2018) Computational Biology and Chemistry, 72, pp. 26-32.
 *
 */
package insilico.coral.models.carcino;


import insilico.coral.CoralAttributes;
import insilico.coral.CoralGraphAttributes;
import insilico.coral.CoralModel;
import insilico.coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 * @author Giovanna
 */
public class CoralFRCarcinogenicity extends CoralModel {

    // Endpoint =  -5.7575532 (± 0.0165842) +    0.1442534 (± 0.0007941) * DCW(2,11)
    public CoralFRCarcinogenicity() {
        super( -5.7575532, 0.1442534, "/carcino_FR_str-cw.txt");
    }

    @Override
    protected double calculateDCW(String smi) throws Exception {
        double DCW = 0;
           
        // Hydrogen suppressed graph (HSG) is used in this model
        // Attributes involve from graph: EC0, NNC,C6, C5
        // Attributes involve from SMILES: S, SS, BOND, NOSP, HALO
        
        // list of tokens from the SMILES in input
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> list_tokens = coral_smi.ProcessSMILES(smi);
        
        // S, SS, BOND, NOSP,HALO.
        ArrayList<String> list_attributes = CoralAttributes.calculateAttributeSks(list_tokens);
        list_attributes.addAll(CoralAttributes.calculateAttributeSSks(list_tokens));
        list_attributes.add(CoralAttributes.calculateAttributeBOND(list_tokens));
        list_attributes.add(CoralAttributes.calculateAttributeNOSPCoral2017(list_tokens));
        list_attributes.add(CoralAttributes.calculateAttributeHALO(list_tokens));
        
        // HSG attributes: EC0, NNC, C5, C6
        
        CoralGraphAttributes coral_graph_attr = new CoralGraphAttributes(smi);
        // C6, C5
        for (int i=6; i>4;i--){
            String AttributeRings = coral_graph_attr.CalculateRings(i);
            list_attributes.add(AttributeRings);
        }        
        // NNC
        list_attributes.addAll(coral_graph_attr.CalculateNNC());
        // EC0
        list_attributes.addAll(coral_graph_attr.CalculateEC(0));
        
/*        for (String str: list_attributes){
            System.out.println(str);
        }
*/        
        ArrayList<Double> weights = calculateSAvalues(list_attributes);
        
        // ADDITIVE scheme for the optimal descriptor
        for (Double d : weights){
            DCW += d;
        }
        
        return DCW;
        
    }
    
}
