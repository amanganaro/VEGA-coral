/**
 * CORAL model for 558 compounds with experimental data for general LOAEL in Rats for oral sub-chronic RDT. 
 * The experimental data is expressed in mg/kg bw (as logarithmic units)
 * The optimal descriptor is based on attributes extracted from SMILES (Sk, SSk, SSSk, APP).
 */
  
package insilico.core.coral.models.loael.general;

import insilico.core.coral.CoralAttributes;
import insilico.core.coral.CoralModel;
import insilico.core.coral.CoralSMILES;

import java.util.ArrayList;

/**
 *
 * @author Giovanna Lavado
 */
public class CoralGeneralLoael extends CoralModel{

    // Endpoint =   2.6738599 (± 0.0057358) +    0.1967440 (± 0.0009514) * DCW(1,30)
    // nuovi valori per l'endpoint
    // 2.359337000 e 0.151302200
    public CoralGeneralLoael(){
        super(2.359337000,0.151302200, "/loael_general_attr-cw_NEW3.txt");
    }



    @Override
    protected double calculateDCW(String smi) throws Exception {
        double DCW = 0;

        // tokens from SMILES string
        CoralSMILES coral_smi = new CoralSMILES();
        ArrayList<String> res = coral_smi.ProcessSMILESwithSquareBracket(smi);

        // Sk,SSk,SSSk on all elements NOT including square brackets
        CoralAttributes.SQUARED_BRACKETS_TOKENS_CALCULATION = true;
        ArrayList<String> attr_list = CoralAttributes.calculateAllLocalAttributes(res);

        // APP: Cl and N, O, S, P
        ArrayList<String> chemicalList = new ArrayList<>();
        chemicalList.add("N");
        chemicalList.add("O");
        chemicalList.add("S");
        chemicalList.add("P");

        // it seems that CORAL uses the ORIGINAL splitting for this purpose, as [N+] and similar token
        // are counted within the couples
        ArrayList<String> tokens = coral_smi.ProcessSMILES(smi);
        ArrayList<String> attrAPPList =
                CoralAttributes.calculateAttributesVectorAPP(tokens, "Cl", chemicalList);

        attr_list.addAll(attrAPPList);

        CoralAttributes.SQUARED_BRACKETS_TOKENS_CALCULATION = false;

        ArrayList<Double> weights = calculateSAvalues(attr_list, false);
        for (Double cw : weights) {
            DCW += cw;
        }


        //////////////////////
//        for (int i=0; i<attr_list.size(); i++) {
//            String at = attr_list.get(i);
//            ArrayList<String> buf = new ArrayList<>();
//            buf.add(at);
//            ArrayList<Double> w = calculateSAvalues(buf);
//
//            System.out.println(attr_list.get(i) + "\t" + (w.size() > 0 ? w.get(0) : "" ) );
//        }
        //////////////////////

        return DCW;
    }

}
