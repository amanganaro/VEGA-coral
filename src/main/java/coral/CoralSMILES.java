package coral;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * 
 * @author Alberto Manganaro
 * @author Giovanna Lavado
 * 
 */
public class CoralSMILES {
    
    private static final String TWO_CHAR_LIST_SOURCE = "/one_element_two_char.txt";
    private final ArrayList<String> TwoCharList;
        

    /**
     * 
     * @throws NullPointerException
     */
    public CoralSMILES() throws Exception {

        // read the list of 2-chars tokens from the local resource
        try {
            URL u = this.getClass().getResource(TWO_CHAR_LIST_SOURCE);
            DataInputStream in = new DataInputStream(u.openStream());
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            
            String line;
            TwoCharList = new ArrayList<>();
            while ((line = br.readLine()) != null){
                TwoCharList.add(line.trim());
            }
            br.close();
        } catch (IOException ex) {
            throw new Exception("Unable to load list of 2-chars");
        }
        
    }
    

    /**
     * Process a SMILES string and return the list of its tokens.
     * 
     * @param SMI
     * @return 
     */
    public ArrayList<String> ProcessSMILES(String SMI) {
        
        ArrayList<String> Tokens = new ArrayList<>();
        
        // HashMap for the substitution of chars:
        // ) to (
        // ] to [
        // . to ^
        HashMap<String, String> map = new HashMap<String, String>();
        map.put(")", "(");
        map.put("]", "[");
        map.put(".", "^"); 
        
        // replace all the needed chars
        for(Map.Entry<String, String> entry : map.entrySet())
            SMI = SMI.replace(entry.getKey(), entry.getValue());
        
        int idx = 0;
        while (idx < SMI.length()) {
        
            // identify current token to be checked to determine if its
            // length is 1, 2 or 3 chars
            String token_to_check = SMI.substring(idx, (idx+3)<SMI.length() ? (idx+3) : SMI.length()  );
            
            // check the current token
            String cur_token = CheckToken(token_to_check, this.TwoCharList);

            // add current token and proceed
            Tokens.add(cur_token);
            idx += cur_token.length();            
        }
        
        return Tokens;
    }
    
    
    /**
     * Check a given token to determine if only a part of it should be 
     * considered a valid token.
     * 
     * @param token
     * @param list
     * @return 
     */
    private String CheckToken(String token, ArrayList<String> list) {
      
        // if token has length = 1, it is automatically a valid token
        if (token.length() == 1)
            return token;
        
        // if token is a cycle id >10 in the form %nn it is a valid token
        if (token.charAt(0)=='%' && token.length()==3) 
            return token;
        
        // check if the first 2 chars of the token should be considered
        // a valid token (using the list of 2-char tokens)
        String two_char = token.substring(0, 2);
        if (TwoCharList.contains(two_char)) 
            return two_char;
        
        // for all remaining cases, the valid token is just the first char
        return "" + token.charAt(0);
    }   
       
}
