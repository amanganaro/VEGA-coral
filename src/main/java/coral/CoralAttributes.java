import java.util.ArrayList;

/**
 * 
 * @author Giovanna Lavado
 * @author Alberto Manganaro
 * 
 */
public class CoralAttributes {

    /**
     * Converts the given string in an integer number according to a conversion in pseudo Base-10. 
     * 
     * @param s a string to convert
     * @return the integer value
     */
    private static int convertBase10(String s){
        int value = 0;
        int max = s.length();
        char[] temp = new char[max];
        int base = 10, pos, dec;
        
        s.getChars(0, max, temp, 0);
        
        for(int k = max; k>=1; k--){
            dec = (int) Math.pow(base, k); // or base^(k-1)?
            pos = max - k;
            value += dec * temp[pos]; // 10^k * ascii_value 
        } 
                
        return value;
    }
    /**
     * Calculate the BOND attribute given the token list.
     * Bond attribute represents the presence or absence of three categories of chemical bonds: 
     * double, triple and stereochemical bonds.      * 
     * Example: BOND11000000
     * 
     * @param tokenList list of tokens.
     * @return the BOND attribute.
     */
    public static String calculateAttributeBOND (ArrayList<String> tokenList){
        String z1 = "BOND";
        String z3 = "0000";
        
        String[][] BONDArray = {
            {"="},
            {"#"},
            {"@", "@@"}}; 
        
        char[] z_2 = {'0','0','0','0'};
    
        int n = BONDArray.length;
        
        for (int i=0; i< n; i++)
            for (String elem : BONDArray[i])
                if (tokenList.contains(elem)){
                    z_2[i]= '1';
                    break;
                }

        String z2 = new String(z_2);
        String attrBOND = z1.concat(z2).concat(z3);
        
        return attrBOND;        
    }
    
    /**
     * Calculate the NOSP attribute given the token list. 
     * NOSP attribute represents the presence or absence of four chemical elements:
     * nitrogen, oxygen, sulphur, phosphorus.
     * 
     * This method is implemented like in Coral version 2017, including some exceptions (Na, Nd, Si, Sm, Pr).
     * Example: NOSP11000000
     * 
     * @param tokenList list of tokens.
     * @return the NOSP attribute.
   */
    public static String calculateAttributeNOSPCoral2017 (ArrayList<String> tokenList){
        String z1 = "NOSP";
        String z3 = "0000";
        
       // chemical elements to search including the exceptions (Na, Nd, Si, Sm, Pr)
        String[][] NOSPArray = {{"N","Na","Nd"},
            {"O"}, 
            {"S", "Si", "Sm"}, 
            {"P","Pr"}};
        
        int n = NOSPArray.length;
        char[] z_2 = {'0','0','0','0'};
        
        for (int i=0; i<n; i++)
            for (String elem : NOSPArray[i])
                if (tokenList.contains(elem)){
                    z_2[i]= '1';
                }               
            
        String z2 = new String(z_2);
        String attrNOSP = z1.concat(z2).concat(z3);
        
        return attrNOSP;
        
    }
    
    /**
     * Calculate the HALO attribute given the token list.
     * HAlO attribute represents the presence or absence of four chemical elements:
     * fluorine, chlorine, bromine and iodine.
     * 
     * This method is implemented like in Coral version 2017, including an exception (In)
     * Example: HALO00100000
     * 
     * @param tokenList list of tokens.
     * @return the HALO attribute
     */
    public static String calculateAttributeHALOCoral2017 (ArrayList<String> tokenList){
        String z1 = "HALO";
        String z3 = "0000";
        
        // chemical elements to search including the exception (In)
        String[][] HALOArray = {{"F"},
            {"Cl"}, 
            {"Br"}, 
            {"I", "In"}}; 
        
        int n = HALOArray.length;
        char[] z_2 = {'0','0','0','0'};
        
        for (int i=0; i<n; i++)
            for (String elem : HALOArray[i])
                if (tokenList.contains(elem)){
                    z_2[i]= '1';
                    break;
                }
        
        String z2 = new String(z_2);
        String attrHALO = z1.concat(z2).concat(z3);
        
        return attrHALO;        
    }
    
    /**
     * Calculate the HARD attribute given the token list. 
     * HARD attribute is an assembling of BOND, NOSP and HALO attributes into a single attribute.
     * 
     * This method is implemented like in Coral version 2017, including some exceptions (Na, Nd, Si, Sm, Pr, In).
     * Example: $00000101000
     * 
     * @param tokenList list of tokens.
     * @return the HARD attribute.
     */
    public static String calculateAttributeHARDCoral2017 (ArrayList<String> tokenList){
        
        // chemical elements to search including the exceptions (Na, Nd, Si, Sm, Pr, In)
        String[][] HARDArray = {{"="}, {"#"}, {"@", "@@"},
            {"N", "Na", "Nd"},  // N 
            {"O"},              // O
            {"S", "Si", "Sm"},  // S
            {"P", "Pr"},        // P
            {"F"}, 
            {"Cl"}, 
            {"Br"}, 
            {"I", "In"}};
        
        int n = HARDArray.length;
        
        String attrHARD = "$00000000000";
        char[] temp = attrHARD.toCharArray();
        
        for (int i=0; i<n; i++){
            for (String element : HARDArray[i])
                if (tokenList.contains(element)){
                    temp[i+1] = '1';
                    break;
                }
        }
        attrHARD = new String(temp);
        
        return attrHARD;
    }

    /**
     * Calculate a list of PAIR attributes given the token list. 
     * PAIR attribute represents the presence of several pairs of chemical elements and/or kind of bonds in a SMILES 
     * (regardless their places in the molecular structure): 
     * fluorine, chlorine, bromine, iodine, nitrogen, oxygen, sulphur, phosphorus and, double and triple bonds.
     * 
     * This method is implemented like in Coral version 2017, including some exceptions (In, Na, Nd, Mn, Sn, Zn, Si, Os, Sm, Pr).
     * Example: ++++N---B2==
     * 
     * @param tokenList list of tokens.
     * @return a list of all PAIR attributes.
     */
    public static ArrayList<String> calculateAttributePAIRCoral2017(ArrayList<String> tokenList){
        String z1 = "++++";
        String attrPAIR = "";
        
        // Chemical elements to serch with their exceptions 
        String[][] PAIRArray = {{"F"},
            {"Cl"}, 
            {"Br"}, 
            {"I", "In"}, 
            {"N", "n", "Na", "Nd", "Mn", "In", "Sn", "Zn"}, 
            {"O"}, 
            {"S", "s", "Si", "Os", "Sm"}, 
            {"P", "Pr"}, 
            {"="}, 
            {"#"}
        };
        
        int n = PAIRArray.length;
        
        ArrayList<String> attrPAIRList = new ArrayList<>();
        int nTokens = tokenList.size();

        // if the number of tokens is greater or equal to two then use two indexes i and j in PAIRArray.  
        // Index i=0,...,n-2 contains the first element to search and j=i+1,...,n-1 contains the second one.
        // Search the pair firstElement, secondElment in the token list and build the attribute
        if (nTokens >=2 ){
            for (int i=0; i<n-1; i++){
                String firstElement = "";
                for (String first : PAIRArray[i])
                if (tokenList.contains(first)){ 
                    for (int j=i+1; j<n; j++){
                        firstElement = PAIRArray[i][0]; // firstElement needs to be set for each iteration of j 
                        String secondElement = "";
                       for (String second : PAIRArray[j])
                        if (tokenList.contains(second)){
                            secondElement = PAIRArray[j][0];
                            // construction of the attribute
                            char[] z_2 = {'-','-','-','-'};
                            char[] z_3 = {'=','=','=','='};
                            
                            // transformations
                            if(firstElement.equals("Cl") && !(secondElement.equals("=") || secondElement.equals("#"))){
                                firstElement = "CL";
                            }
    
                            if (firstElement.equals("="))
                                firstElement = "B2";
                            if (secondElement.equals("="))
                                secondElement = "B2";
                            
                            if (firstElement.equals("#"))
                                firstElement = "B3";
                            if (secondElement.equals("#"))
                                secondElement = "B3";
                            
                            // firstElement and secondElement in z_2 and z_3
                            firstElement.getChars(0, firstElement.length(), z_2, 0);
                            secondElement.getChars(0, secondElement.length(), z_3, 0);
                            String z2 = new String(z_2);
                            String z3 = new String(z_3);
                            
                            // concatenation of z1, z2 and z3
                            attrPAIR = z1.concat(z2).concat(z3);
                            attrPAIRList.add(attrPAIR);
                            break;
                        }    
                    }
                    break;    
                }
            }    
        }
        
        return attrPAIRList;
    }
    
    /**
     * Calculate SSSk attributes given the token list.
     * Example: Br..n...F...
     * 
     * @param tokenList list of tokens.
     * @return a list of all SSSk attributes.
     */
    public static ArrayList<String> calculateAttributeSSSks(ArrayList<String> tokenList){
        
        ArrayList<String> attrSSSkList = new ArrayList<>();
        
        String attrSSSk = "";
        
        char[] z_1 = {'.','.','.','.'};
        char[] z_2 = {'.','.','.','.'};
        char[] z_3 = {'.','.','.','.'};
        
        int idx = 0;
        int nTokens = tokenList.size();
        
        if (nTokens >=3 ){
            while (idx < nTokens-2){

                String token1 = tokenList.get(idx);
                String token2 = tokenList.get(idx+1);
                String token3 = tokenList.get(idx+2);
                String temp = null;
                
                // initialization z_1, z_2, z_3
                for (int j=0; j<4;j++) {
                    z_1[j] = '.';
                    z_2[j] = '.';
                    z_3[j] = '.';
                }
             
                token1.getChars(0, token1.length(), z_1, 0);
                token2.getChars(0, token2.length(), z_2, 0);
                token3.getChars(0, token3.length(), z_3, 0);
                
                String z1 = new String(z_1);
                String z2 = new String(z_2);
                String z3 = new String(z_3);

                // decreasing order between z1 and z3
                if (convertBase10(z1) < convertBase10(z3)){
                    temp = z1;
                    z1 =  z3;
                    z3 = temp;
                }
                
                // concatenation z1,z2 and z3
                attrSSSk = z1.concat(z2).concat(z3);
                
                attrSSSkList.add(attrSSSk);
                
                idx ++;
            }
        
        }

        return attrSSSkList;
    }
    
    
    /**
     * Calculate SSk attributes given the token list.
     * Example: Br..F.......
     * 
     * @param tokenList a list of tokens.
     * @return a list of all SSk attributes.
     */
    public static ArrayList<String> calculateAttributeSSks(ArrayList<String> tokenList){
        String z3 = "....";
        String attrSSk = "";
        
        ArrayList<String> attrSSkList = new ArrayList<>();
    
        int idx = 0;
        int nTokens = tokenList.size();
        
        if (nTokens >=2 ){
            while (idx < nTokens-1){ 

                String token1 = tokenList.get(idx);
                String token2 = tokenList.get(idx+1);
                String temp = null;
                
                char[] z_1 = {'.','.','.','.'};
                char[] z_2 = {'.','.','.','.'};
                    
                token1.getChars(0, token1.length(), z_1, 0);
                token2.getChars(0, token2.length(), z_2, 0);
                
                String z1 = new String(z_1);
                String z2 = new String(z_2);
                
                // ordering between z1 and z2, decreasing order
                if (convertBase10(z1) < convertBase10(z2)){
                    temp = z1;
                    z1 = z2;
                    z2 = temp;
                }
                
                // concatenation z1,z2 and z3
                attrSSk = z1.concat(z2).concat(z3);
                
                attrSSkList.add(attrSSk);
                
                idx ++;
            }
        
        }
            
        return attrSSkList;
    }
    
    
    /**
     * Calculate Sk attributes given the token list.  
     * Example: Cl..........
     * 
     * @param tokenList a list of tokens.
     * @return a list of all Sk attributes.
     */
    public static ArrayList<String> calculateAttributeSks(ArrayList<String> tokenList){
        String z2 = "....";
        String z3 = "....";
        
        String attrSk = "";
        
        ArrayList<String> attrSkList = new ArrayList<>();
        
        // each token is copied in z 
        for (String token : tokenList){
            char[] z = {'.','.','.','.'};
            token.getChars(0, token.length(), z, 0);
            
            String z1 = new String(z);
            
            // the atttribute is the concatenation of z1, z2 and z3
            attrSk = z1.concat(z2).concat(z3);
            
            attrSkList.add(attrSk);
            
        }
              
        return attrSkList;
    }
    
    
    /**
     * Calculate all the Sk, SSk, SSSk attributes for the given tokens list.
     * 
     * @param tokenList a list of tokens.
     * @return a list of all local attributes.
     */
    public static ArrayList<String> calculateAllLocalAttributes(ArrayList<String> tokenList){
        
        ArrayList<String> localAttrList = calculateAttributeSks(tokenList);
        localAttrList.addAll(calculateAttributeSSks(tokenList));
        localAttrList.addAll(calculateAttributeSSSks(tokenList));
        
        return localAttrList;
    }
    
    /**
     * Calculate the HALO attribute given the token list.
     * HAlO attribute represents the presence or absence of four chemical elements:
     * fluorine, chlorine, bromine and iodine.
     * Example: HALO00100000
     * 
     * @param tokenList list of tokens.
     * @return the HALO attribute
     */    
    public static String calculateAttributeHALO (ArrayList<String> tokenList){
        String z1 = "HALO";
        String z3 = "0000";
        
        // chemical elements to search: F, Cl, Br, I
        String[] HALOArray = {"F", "Cl", "Br", "I"};
        
        char[] z_2 = {'0','0','0','0'};
        int i = 0;
        for (String elem : HALOArray){
            if (tokenList.contains(elem))
                z_2[i]= '1';
            i++;
        }
        
        String z2 = new String(z_2);
        String attrHALO = z1.concat(z2).concat(z3);
        
        return attrHALO;
    }
    
    /**
     * Calculate the NOSP attribute given the token list. 
     * NOSP attribute represents the presence or absence of four chemical elements:
     * nitrogen, oxygen, sulphur, phosphorus.
     * Example: NOSP11000000
     * 
     * @param tokenList list of tokens.
     * @return the NOSP attribute.
 
     */
    public static String calculateAttributeNOSP (ArrayList<String> tokenList){
        String z1 = "NOSP";
        String z3 = "0000";
        
        // chemical elements to search
        String[] NOSPArray = {"N", "O", "S", "P"};
        
        char[] z_2 = {'0','0','0','0'};
        int i = 0;
        for (String elem : NOSPArray){
            if (tokenList.contains(elem))
                z_2[i]= '1';
            i++;
        }
        
        String z2 = new String(z_2);
        String attrNOSP = z1.concat(z2).concat(z3);
        
        return attrNOSP;
    }
    
    /**
     * Calculate the HARD attribute given the token list. 
     * HARD attribute is an assembling of BOND, NOSP and HALO attributes into a single attribute.
     * Example: $00000101000
     * 
     * @param tokenList list of tokens.
     * @return the HARD attribute.
     */
    public static String calculateAttributeHARD (ArrayList<String> tokenList){
        
        String[][] HARDArray = {{"="}, {"#"}, {"@", "@@"}, // BOND
            {"N"},  // N 
            {"O"},  // O
            {"S"},  // S
            {"P"},  // P
            {"F"}, // HALO -->
            {"Cl"}, 
            {"Br"}, 
            {"I"}};
        
        int n = HARDArray.length;
        
        String attrHARD = "$00000000000";
        char[] temp = attrHARD.toCharArray();
        
        for (int i=0; i<n; i++){
            for (String element : HARDArray[i])
                if (tokenList.contains(element)){
                    temp[i+1] = '1';
                    break;
                }
        }
        attrHARD = new String(temp);
        
        return attrHARD;
    }
    
       /**
     * Calculate a list of PAIR attributes given the token list. 
     * PAIR attribute represents the presence of several pairs of chemical elements and/or kind of bonds in a SMILES 
     * (regardless their places in the molecular structure): 
     * fluorine, chlorine, bromine, iodine, nitrogen, oxygen, sulphur, phosphorus and, double and triple bonds.
     * Example: ++++N---O===
     * 
     * @param tokenList list of tokens.
     * @return a list of all PAIR attributes.
     */
    public static ArrayList<String> calculateAttributePAIR(ArrayList<String> tokenList){
        String z1 = "++++";
        String attrPAIR = "";
        
        // chemical elements and/or kind of bonds to search        
        String[] PAIRArray = {"F", "Cl", "Br", "I", "N", "O", "S", "P", "=", "#"};
        int n = PAIRArray.length;

        ArrayList<String> attrPAIRList = new ArrayList<>();
        int nTokens = tokenList.size();
       
        // if the number of tokens is greater or equal to two then, use two indexes i and j in PAIRArray.  
        // Index i=0,...,n-2 contains the first element to search and j=i+1,...,n-1 contains the second one.
        // Search the pair firstElement, secondElment in the token list and build the attribute
        if (nTokens >=2 ){
            for (int i=0; i<n-1; i++){
                String firstElement = "";
                if (tokenList.contains(PAIRArray[i])){
                    firstElement = PAIRArray[i];
                    for (int j=i+1; j<n; j++){
                        String secondElement = "";
                        if (tokenList.contains(PAIRArray[j])){
                            secondElement = PAIRArray[j];

                            // building the attribute
                            char[] z_2 = {'-','-','-','-'};
                            char[] z_3 = {'=','=','=','='};
                            
                            // transformations
                            if (firstElement.equals("="))
                                firstElement = "B2";
                            if (secondElement.equals("="))
                                secondElement = "B2";
                            
                            if (firstElement.equals("#"))
                                firstElement = "B3";
                            if (secondElement.equals("#"))
                                secondElement = "B3";
                            
                            // first and second element in z_2 e z_3, respectively
                            firstElement.getChars(0, firstElement.length(), z_2, 0);
                            secondElement.getChars(0, secondElement.length(), z_3, 0);
                            String z2 = new String(z_2);
                            String z3 = new String(z_3);
                            
                            // concatenation of z1, z2 and z3
                            attrPAIR = z1.concat(z2).concat(z3);
                            attrPAIRList.add(attrPAIR);
                        }    
                    }  
                }
            }    
        }
          
        return attrPAIRList;
    }

    public static String calculateAttributesPairProportions(ArrayList<String> tokenList, String ChemicalOne, String ChemicalTwo){
        String z3 = "....";
         
        int n = countOccurrencesOfChemical(tokenList, ChemicalOne);
        if (n > 9) n = 9;
                      
        int m = countOccurrencesOfChemical(tokenList, ChemicalTwo);
        if (m > 9) m = 9;
        
        String attr = "(" + ChemicalOne + ChemicalTwo + ")";
        if (ChemicalTwo.equals("Cl")) attr = attr + n;
        else attr = attr + "." + n;
        
        attr = attr + "." + m + z3;
                
        return attr;
       
    }
    
    private static int countOccurrencesOfChemical(ArrayList<String> tokenList, String Chemical){
        int retval = 0;

        String[] myArr = new String[tokenList.size()];
        myArr = tokenList.toArray(myArr); 
        int nTokens = myArr.length;
        
        for (int i=0 ; i<nTokens; i++)
            if (myArr[i].equals(Chemical)){
                retval++;
        }
                
        return retval;
    }    
    
}
