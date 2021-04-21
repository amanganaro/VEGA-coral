package insilico.coral;

import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InvalidMoleculeException;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;

import java.util.ArrayList;

/**
 * Calculation of the CORAL graph attributes.
 * 
 * @author Giovanna Lavado
 * @author Alberto Manganaro
 */
public class CoralGraphAttributes {
    
    final InsilicoMolecule CurMol;
    final double[][] AdjMatrix;
    final IRingSet CurRings;
    
    
    /**
     * Constructor. Convert the given molecule from a SMILES string to an
     * internal format where the molecular graph is available. Calculate
     * the adjacency connection matrix. Calculate the set of all rings found
     * in the molecular structure.
     *
     * All the computations are based on a H-depleted molecule.
     * 
     * @param SMILES
     * @throws GenericFailureException 
     * @throws InvalidMoleculeException
     * @throws CDKException
     */
    public CoralGraphAttributes(String SMILES) throws GenericFailureException,
            InvalidMoleculeException, CDKException {
        
        // convert the given SMILES string to a InsilicoMolecule object
        // which contains the molecule as a graph, together with other 
        // chemical information
        CurMol = SmilesMolecule.Convert(SMILES);
        
        // if there is a parsing error (the SMILES string is not correct)
        // throws an exception as no further elaborations can be done
        if (!CurMol.IsValid())
            throw new GenericFailureException(CurMol.GetErrors().GetMessages());
        
        // get the adjacency connection matrix, with some modification:
        // 1. aromatic bonds are set to 1 (originally coded as 1.5)
        // 2. diagonal entries are set to 0 
        AdjMatrix = CurMol.GetMatrixConnectionAugmented();
        for (int i=0; i<AdjMatrix.length; i++)
            for (int j=0; j<AdjMatrix[0].length; j++) {
                if (i==j) 
                    AdjMatrix[i][j] = 0;
                if (AdjMatrix[i][j] == 1.5)
                    AdjMatrix[i][j] = 1;
            }
        
        // retrieve the full set of all rings available in the current molecule
        AllRingsFinder ARF = new AllRingsFinder();
        CurRings =  ARF.findAllRings(CurMol.GetStructure());        
    }
    
    /**
     * Calculate the ring attributes for a given ring size.
     * 
     * @param RingSize number of atoms for the requested ring
     * @return the retrieved attribute
     */
    public String CalculateRings(int RingSize) {
        
        // total number of rings of the current size
        int RingCount = 0;
        
        // flag, if at least one ring is an aromatic rings
        boolean hasAromaticRings = false;

        // flag, if at least one ring has hetero-atoms
        boolean hasRingsWithHeteroAtoms = false;
        
        // iteration for all found rings, check their size (= number of atoms)
        for (IAtomContainer CurRing : CurRings.atomContainers()) {
            int nAtoms = CurRing.getAtomCount();
            
            if (nAtoms == RingSize) {
                
                // found a ring of the requested size. 
                RingCount++;
                
                // check for additional information for the attribute
                boolean isCurRingAromatic = true;
                for (IAtom CurAtom : CurRing.atoms()) {
                    if (!CurAtom.getSymbol().equalsIgnoreCase("C"))
                        hasRingsWithHeteroAtoms = true;
                    if (!CurAtom.getFlag(CDKConstants.ISAROMATIC))
                        isCurRingAromatic = false;
                } 
                if (isCurRingAromatic)
                    hasAromaticRings = true;

            }
                        
        }

        // build and return the attribute 
        String z_1 = "C" + RingSize + "..";
        String z_2 = "." + (hasAromaticRings ? "A" : ".") +
                (hasRingsWithHeteroAtoms ? "H" : ".") + ".";
        String z_3 = String.valueOf(RingCount) + (RingCount > 9 ? ".." : "...");
        String CurAttribute = z_1 + z_2 + z_3;
        
        return CurAttribute;                
    }


    /**
     * Calculate the ring attributes for a given ring size (only hetero-atoms are considered).
     *
     * @param RingSize number of atoms for the requested ring
     * @return the retrieved attribute
     */

    public String CalculateRingsWithHeteroAtoms(int RingSize) {

        // total number of rings of the current size
        int RingCount = 0;

        // flag, if at least one ring has hetero-atoms
        boolean hasRingsWithHeteroAtoms = false;

        // iteration for all found rings, check their size (= number of atoms)
        for (IAtomContainer CurRing : CurRings.atomContainers()) {
            int nAtoms = CurRing.getAtomCount();

            if (nAtoms == RingSize) {

                // found a ring of the requested size.
                RingCount++;

                for (IAtom CurAtom : CurRing.atoms()) {
                    if (!CurAtom.getSymbol().equalsIgnoreCase("C"))
                        hasRingsWithHeteroAtoms = true;
                }

            }

        }

        // build and return the attribute
        String z_1 = "C" + RingSize + "..";
        String z_2 = ".." +
                (hasRingsWithHeteroAtoms ? "H" : ".") + ".";
        String z_3 = String.valueOf(RingCount) + (RingCount > 9 ? ".." : "...");
        String CurAttribute = z_1 + z_2 + z_3;

        return CurAttribute;
    }
    
    /**
     * Calculate the NNC (Nearest Neighbours Code) attributes.
     * 
     * @return list of the retrieved attributes
     * @throws InvalidMoleculeException 
     */
    public ArrayList<String> CalculateNNC() throws InvalidMoleculeException {
    
        // number of atoms, also the size of the adjacency matrix (square matrix)
        int nAt = AdjMatrix.length;
        
        // iteration on all atoms, to calculate the NNC for each of them
        ArrayList<String> Attributes = new ArrayList<>();
        for (int atom=0; atom<nAt; atom++) {
            
            String CurAtomName = CurMol.GetStructure().getAtom(atom).getSymbol();
            
            // number of connected atoms to the current atom
            int nNeighbours = 0;
            
            // number of connected atoms that are Carbon atoms
            int nCarbonNeighbours = 0;

            // number of connected atoms that are non-Carbon atoms
            int nNonCarbonNeighbours = 0;
            
            // iteration on all atoms, to check if they are connected to
            // the current atom
            for (int i=0; i<nAt; i++) 
                if (AdjMatrix[atom][i] > 0) {
                    
                    // increment number of neighbours
                    nNeighbours++;
                    
                    // retrieve name (symbol) of connected atom
                    String ConnectedAtomName = CurMol.GetStructure().getAtom(i).getSymbol();
                    
                    // check if the connected atom is a carbon
                    if (ConnectedAtomName.equalsIgnoreCase("C"))
                        nCarbonNeighbours++;
                    else
                        nNonCarbonNeighbours++;
                }
            
            // GIOVANNA: costruisci qui il nome dell'attributo, non ricordo
            // bene il formato, lascio un codice simile a quello che dovrai
            // fare tu
            // CurAtomName occupa i successivi 4 caratteri ini ....
            // uguale dopo
            char[] z_2 = {'.','.','.','.'};
            char[] z_3 = {'.','.','.','.'};
            
            CurAtomName.getChars(0, CurAtomName.length(), z_2, 0);
            String z2 = new String(z_2);
            
            String NN_code = String.valueOf(nNeighbours) +
                    String.valueOf(nCarbonNeighbours) + String.valueOf(nNonCarbonNeighbours);
            
            NN_code.getChars(0, NN_code.length(), z_3, 0);
            String z3 = new String(z_3);
            
            //String CurNNC = "NNC-" + CurAtomName + ".." + nNeighbours +
            //        nCarbonNeighbours + nNonCarbonNeighbours;
            String CurNNC = "NNC-" + z2 + z3;
            Attributes.add(CurNNC);
            
        }
        
        return Attributes;
    }
    
    
    
    /**
     * Calculate the EC (Extended connectivity) attributes.
     * 
     * @param MaxLag max lag for the EC calculation
     * @return list of the retrieved attributes
     * @throws InvalidMoleculeException 
     */
    public ArrayList<String> CalculateEC(int MaxLag) throws InvalidMoleculeException {
        
        // number of atoms, also the size of the adjacency matrix (square matrix)
        int nAt = AdjMatrix.length;
        
        // prepare the data structure for the results
        int[][] ECList = new int[MaxLag+1][nAt];
        
        // calculate the EC0
        for (int atom=0; atom<nAt; atom++) {
            
            // reset value for current EC0
            ECList[0][atom] = 0;
            
            // check all connection for current atom in the adj matrix
            // for each connected atom, EC0 is incremented by 1
            for (int i=0; i<nAt; i++) 
                if (AdjMatrix[atom][i] > 0)
                    ECList[0][atom]++;
            
        }
        
        // calculate iteratively the next EC lags
        for (int lag=1; lag<=MaxLag; lag++) {
            
            for (int atom=0; atom<nAt; atom++) {

                // reset value for current EC lag
                ECList[lag][atom] = 0;

                // check all connection for current atom in the adj matrix
                // for each connected atom, current EC is incremented by the
                // value of the previous EC lag of the connected atom
                for (int i=0; i<nAt; i++) 
                    if (AdjMatrix[atom][i] > 0)
                        ECList[lag][atom] +=  ECList[lag-1][i];

            }
            
        }
        
        // build the attributes
        ArrayList<String> Attributes = new ArrayList<>();
        for (int lag=0; lag<=MaxLag; lag++) {
            for (int atom=0; atom<nAt; atom++) {
                
                String CurECName = "EC" + lag;
                String CurAtomName = CurMol.GetStructure().getAtom(atom).getSymbol();
                String CurValue = Integer.toString(ECList[lag][atom]);
                
                // GIOVANNA: costruisci qui il nome dell'EC corrente, non ricordo
                // bene il formato, lascio un codice simile a quello che dovrai
                // fare tu
		// CurAtomName occupa i sucessivi 4 caratteri inizializzati ....
                // CurVAlue anche
                char[] z_2 = {'.','.','.','.'};
                char[] z_3 = {'.','.','.','.'};

                CurAtomName.getChars(0, CurAtomName.length(), z_2, 0);
                CurValue.getChars(0, CurValue.length(), z_3, 0);
                
                String z2 = new String(z_2);
                String z3 = new String(z_3);
                
                String CurEC = CurECName + "-" + z2 + z3;
                Attributes.add(CurEC); 
                
            }
        }
        
        return Attributes;
    }
    
}
