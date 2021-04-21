package insilico.coral;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;

/**
 *
 * @author Giovanna Lavado
 * @author Alberto Manganaro
 * 
 */
public abstract class CoralModel {

    // activity = c0 (+/- e0) + c1 (+/- e1)* dcw(t^*,n^*)
    
    protected final double intercept; // the intercept is c0 in coral
    protected final double slope; // the slope is c1 in coral
    protected final String SaSource;
    protected ArrayList<StructuralAttribute> SAsList; // structural attributes: attr:cw
    

    public CoralModel(double intercept, double slope, String SaSource) {
        this.intercept = intercept;
        this.slope = slope;
        this.SaSource = SaSource;
        this.SAsList = null;
    }
  
    
    /**
     * Read the SAs list and weights from text file provided as local resource
     * (URL given in the constructor).
     * 
     * @throws Exception
     */
    protected void InitSAs() throws Exception {
        // read SA list and weights from the given resource
        URL u = this.getClass().getResource(SaSource);
        this.SAsList = readStructuralAttributeFromFile(u.openStream(),1);        
    }
    
    
    /**
     * Generic function for DCW calculation, to be implemented in the subclass
     * for each specific model.
     * 
     * @param smi
     * @return 
     * @throws Exception
     */
    protected abstract double calculateDCW(String smi) throws Exception;
    
    
    /**
     * Retrieve the list of weights for all the given SAs, from the list of
     * current SAs weights loaded in the constructor. Values are provided as a
     * list so that different schemes (additive etc.) can be applied in the 
     * model.
     * 
     * @param SAs
     * @return 
     */
    protected ArrayList<Double> calculateSAvalues(ArrayList<String> SAs) {
    
        ArrayList<Double> Values = new ArrayList<>();
        
        for (String attrInSmiles : SAs)
            for (StructuralAttribute str_attr : this.SAsList) 
                if (attrInSmiles.compareTo(str_attr.getAttribute())==0)
                    Values.add(str_attr.getCorrWeight());
     
        return Values;
    }
    
    
    /**
     * Return the prediction for a given SMILES structure, based on the 
     * parameters and weights provided to the constructor. Calculation of DCW
     * for the molecule relies on the class implemented in the specific
     * model subclass.
     * 
     * @param smi
     * @return 
     */
    public double Predict(String smi) throws Exception {
        
        // check if weight have been loaded
        if (SAsList == null)
            InitSAs();
        
        // calculate prediction
        return this.intercept + this.slope*calculateDCW(smi);
    }
    
        
    /**
     * Read the list of SA and weight from the given resource, starting from
     * the given line number. 
     * 
     * @param source
     * @param startLine
     * @return
     * @throws Exception
     */
    private static ArrayList<StructuralAttribute> readStructuralAttributeFromFile(InputStream source,
                                                                                  int startLine) throws Exception {
        
        ArrayList<StructuralAttribute> SAsList = new ArrayList<>();
        
        int currentLine = 0;
        
         try {
             
            DataInputStream in = new DataInputStream(source);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
             
             // skip firts startLine-1 rows
             while(currentLine < startLine-1){
                 if (br.readLine()==null) throw new IOException("File too small: reading row " + startLine);
                 currentLine ++;
             }
             
             // read the number of attributes that should be located on the
             // first valid row of the file
             String line = br.readLine();
             int nAttr = Integer.parseInt(line);
             currentLine++;
             
            // read all the following nAttr lines for the attributes
            while (currentLine < startLine + nAttr){
                if ((line=br.readLine())==null) throw new IOException("File too small: reading row " + currentLine);

                // source file should be tab-separed plain text
                String[] data = line.split("\t");
                
                // source file should have the attribute in the first column
                // and its weight in the second column
                StructuralAttribute attr = createStructuralAttribute(data);
                
                if (attr != null) 
                    SAsList.add(attr);
                
                currentLine++;
            }
             
            br.close();
          
            if (SAsList.size() != nAttr)
                throw new Exception("number of attributes not matching with the number declared in the file");
                         
            return SAsList;             
        
        } catch (IOException | NumberFormatException e) {
            throw new Exception("Unable to read from SAs source - " + e.getLocalizedMessage());
        } 
    }    


    /**
     * Create a single SA with the given pair of attribute-value.
     * 
     * @param data
     * @return
     * @throws Exception
     */
    private static StructuralAttribute createStructuralAttribute(String[] data) throws Exception {
        String attribute = data[0];
        try {
            double corrWeight = Double.parseDouble(data[1]);
            return new StructuralAttribute(attribute, corrWeight);
        } catch (NumberFormatException ex){
            throw new Exception("unable to convert numerical value for input: " + data[1]);
        }
    }    
    
}