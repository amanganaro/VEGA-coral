package insilico.core.coral;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.core.molecule.matrix.AdjacencyMatrix;
import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;

import static insilico.core.molecule.tools.Manipulator.CountImplicitHydrogens;

/**
 * GAO descriptors.
 * 
 * @author Alberto Manganaro
 */
public class GAODescriptor {

    public final static short H_TYPE_ALL = 1; // add H the standard way
    public final static short H_TYPE_HETERO_TERMINAL = 2; // add H only on hetero and terminal
    public final static short H_TYPE_AROMATIC = 3; // add H only on aromatic atoms
    
    private final HashMap<String, ArrayList<String>> OrbitalForAtomType;
    
    private int[][] GAOMatrix;
    private String[] GAOMatrixHeaders;
    private short HAddType;

    
    
    /**
     * Constructor, init the list of orbitals available for known atom types.
     * Orbitals list taken from Torpov et al. "OCWLGI Descriptors: Theory and 
     * Praxis", 2013.
     * 
     * @throws IOException
     * @throws insilico.core.exception.InitFailureException
     */
    public GAODescriptor() throws IOException, InitFailureException {
        
        GAOMatrix = null;
        GAOMatrixHeaders = null;
        HAddType = H_TYPE_ALL;

        // Init list of atomic orbitals for known atom types
        // Note: Bk should be 1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10  4f14 5s2 5p6 5d10  5f9 6s2 6p6 7s2.
        // But we have kept the original file
        URL u = getClass().getResource("/AOs.txt");
        this.OrbitalForAtomType = readAOListFromFile(u.openStream());
        
    }
    
    
    /**
     * Overload of Calculate() method, to be used if only the SMILES is 
     * available (slower, it has to parse and convert the molecule).
     * 
     * @param SMILES 
     * @throws InvalidMoleculeException 
     * @throws insilico.core.exception.GenericFailureException 
     */
    public void Calculate(String SMILES) throws
            InvalidMoleculeException, GenericFailureException {
        Calculate(SmilesMolecule.Convert(SMILES));
    }
    
    
    /**
     * Calculate the GAO matrix.
     * 
     * @param Mol 
     * @throws InvalidMoleculeException
     * @throws GenericFailureException 
     */
    public void Calculate(InsilicoMolecule Mol) throws 
            InvalidMoleculeException, GenericFailureException {

        // Check if the molecule is ok
        if (!Mol.IsValid())
            throw new InvalidMoleculeException("Invalid molecule");
        
        // GAO descriptors work on H-filled structures, add hydrogen
        // Use a custom function to be compliant with CORAL
        IAtomContainer curMol = AddHydrogensCustom(Mol.GetStructure());
        
        // Calculate adjacency matrix
        int[][] AdjMat = AdjacencyMatrix.getMatrix(curMol);
        
        
        //// Calculate GAO matrix from adjacency matrix
        
        // Build list of orbitals for each molecule's atom, following the 
        // atom order/numbering of the molecule
        
        int nAtoms = curMol.getAtomCount();
        ArrayList<String>[] OrbitalsForAtom = new ArrayList[nAtoms];
        int nTotOrbitals = 0;
        
        for (int atIndex = 0; atIndex < nAtoms; atIndex++) {
            
            // Get atom type of current atom
            IAtom curAtom = curMol.getAtom(atIndex);
            String curAtomType = curAtom.getSymbol();
            
            // Retrieve orbitals for the atom and store it in the array
            if (!OrbitalForAtomType.containsKey(curAtomType))
                throw new GenericFailureException("Orbitals list not found for atom type: " + curAtomType);
            ArrayList<String> curOrbitals = OrbitalForAtomType.get(curAtomType);
            OrbitalsForAtom[atIndex] = curOrbitals;
            nTotOrbitals += curOrbitals.size();
            
        }
        
        // Build matrix and headers (with orbital names)
        
        GAOMatrix = new int[nTotOrbitals][nTotOrbitals];
        GAOMatrixHeaders = new String[nTotOrbitals];
        int idxGAOMatrix = 0;
        
        for (int atIndex = 0; atIndex < nAtoms; atIndex++) {
            ArrayList<String> curOrbitals = OrbitalsForAtom[atIndex];
            for (int orbIndex = 0; orbIndex < curOrbitals.size(); orbIndex++) {
                
                // Set string for headers (orbital name)
                GAOMatrixHeaders[idxGAOMatrix] = curOrbitals.get(orbIndex);
                
                // Run through all other atoms and fill the matrix
                int idxGAOMatrix2 = 0;
                for (int atIndex2 = 0; atIndex2 < nAtoms; atIndex2++) {
                    
                    ArrayList<String> curOrbitals2 = OrbitalsForAtom[atIndex2];
                    for (int orbIndex2 = 0; orbIndex2 < curOrbitals2.size(); orbIndex2++) {
                        
                        // if on diagonal of adjacency matrix, always zero
                        if (atIndex == atIndex2)
                            GAOMatrix[idxGAOMatrix][idxGAOMatrix2] = 0;
                        
                        // if adjacency matrix = 0, also GAO = 0
                        else if (AdjMat[atIndex][atIndex2] == 0) 
                            GAOMatrix[idxGAOMatrix][idxGAOMatrix2] = 0;
                        
                        // else if adjacency matrix = 1, also GAO = 1
                        else
                            GAOMatrix[idxGAOMatrix][idxGAOMatrix2] = 1;
                        
                        idxGAOMatrix2++;
                    }
                    
                }
                
                idxGAOMatrix++;
                                
            }
        }
        
    }
 
    
    /**
     * 
     * @param type 
     */
    public void SetHAddType(short type) {
        this.HAddType = type;
    }
    
    
    /**
     * @return the GAOMatrix
     */
    public int[][] getGAOMatrix() {
        return GAOMatrix;
    }

    
    /**
     * @return the GAOMatrixHeaders
     */
    public String[] getGAOMatrixHeaders() {
        return GAOMatrixHeaders;
    }
    
    
    /**
     * Utility to print (to a output stream like System.out) the calculated
     * GAO matrix.
     * 
     * @param out 
     */
    public void PrintGAOMatrix(PrintStream out) {
        if (GAOMatrix == null)
            return;
        
        // First row, headers
        out.print("/");
        for (int i=0; i<GAOMatrix.length; i++) 
            out.print("\t" + GAOMatrixHeaders[i]);
        out.println();
        
        // Matrix
        for (int i=0; i<GAOMatrix.length; i++) {
            out.print(GAOMatrixHeaders[i]);
            for (int j=0; j<GAOMatrix.length; j++) {
                out.print("\t" + GAOMatrix[i][j]);
            }     
            out.println();
        }
    }
    
    
    /**
     * Custom version of AddHydrogens(). It is needed because apparently in 
     * CORAL the molecule is NOT H-filled, but H are added only on some
     * kind of atoms. The calculation should be compliant with CORAL.
     * Original class: insilico.core.molecule.tools.Manipulator
     * 
     * @param mol
     * @return
     * @throws GenericFailureException 
     */
    private IAtomContainer AddHydrogensCustom(IAtomContainer mol)
        throws GenericFailureException {

        IAtomContainer NewMol=null;

        try {
            NewMol= mol.clone();
        } catch (CloneNotSupportedException e) {
            throw new GenericFailureException("Unable to clone molecule");
        }

        int OldCount = mol.getAtomCount();

        for (int i=0;i<OldCount;i++) {
            IAtom atom = NewMol.getAtom(i);
            
            // check if hetero or terminal atom
            if (HAddType == H_TYPE_HETERO_TERMINAL) {
                boolean hetero = false;
                boolean terminal = false;
                if (!atom.getSymbol().equalsIgnoreCase("C"))
                    hetero = true;
                if (NewMol.getConnectedAtomsCount(atom) == 1)
                    terminal = true;
                if ( (!hetero) && (!terminal) ) 
                    continue;
            }

            // check if aromatic
            if (HAddType == H_TYPE_AROMATIC) {
                if (!atom.getFlag(CDKConstants.ISAROMATIC))
                    continue;
            }
            
            
            int HCount = CountImplicitHydrogens(atom);

            atom.setImplicitHydrogenCount(0);

            for (int j=1;j<=HCount;j++) {
                Atom hydrogen = new Atom("H");
                hydrogen.setAtomicNumber(1);
                NewMol.addAtom(hydrogen);
                Bond newBond = new Bond(atom, hydrogen, IBond.Order.SINGLE);
                NewMol.addBond(newBond);
            }
            
        }

        return NewMol;
    } 
    
    private static HashMap<String, ArrayList<String>> readAOListFromFile(InputStream source) throws InitFailureException, IOException {
        HashMap<String, ArrayList<String>> OrbitalForAtomType = new HashMap<>();
        ArrayList<String> orbitals = new ArrayList<>();
        
        DataInputStream in = new DataInputStream(source);
        BufferedReader br = new BufferedReader(new InputStreamReader(in));
        
        // read the first line from the text file
        String line = br.readLine();
        
        // loop until all lines are read
            while (line != null){
                // source file should be space-separed plain text
                String[] data = line.split(" ");
                
                // source file should have symbol and a list of orbitals
                orbitals = new ArrayList<>();
                for (int i = 1; i < data.length; i++) {
                    orbitals.add(data[i]);
                }
                OrbitalForAtomType.put(data[0], orbitals);    
                line = br.readLine();
            }
        
        br.close();
        
        if (OrbitalForAtomType.isEmpty())
            throw new InitFailureException("AO list is empty");
        
        return OrbitalForAtomType;
    }
    
    /**
     * Calculate the EC (Extended connectivity) attributes of zero order (vertex degree).
     * 
     * @return list of the retrieved attributes
     */
    public ArrayList<String> CalculateEC0(){
        ArrayList<String> GAODescriptorsList = new ArrayList<>();
        
        for (int i=0; i<GAOMatrix.length; i++) {
            
            // Calculate EC0 (sum of the row) for the current orbital
            int EC = 0;
            for (int j=0; j<GAOMatrix.length; j++)
                EC += GAOMatrix[i][j];
            
            // Build descriptor name
            String DescName = "EC0-" + GAOMatrixHeaders[i] + " ";
            DescName += Integer.toString(EC);
            if (EC<10)
                DescName += "...";
            else if (EC<100)
                DescName += "..";
            else if (EC<1000)
                DescName += ".";
            
            GAODescriptorsList.add(DescName);
        }
        
        return GAODescriptorsList;
    
    }
    
    /**
     * Calculate the EC (Extended connectivity) attributes.
     * 
     * @param MaxLag max lag for the EC calculation
     * @return list of the retrieved attributes
     */
    public ArrayList<String> CalculateEC(int MaxLag){
        ArrayList<String> GAODescriptorsList = new ArrayList<>();
                
        return GAODescriptorsList;
    }
    
    /**
     * Calculate the NNC (Nearest Neighbours Code) attributes.
     * 
     * @return list of the retrieved attributes
     */
    public ArrayList<String> CalculateNNC(){
        ArrayList<String> GAODescriptorsList = new ArrayList<>();
                
        return GAODescriptorsList;
    }
}
