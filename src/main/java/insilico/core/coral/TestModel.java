package insilico.core.coral;

import insilico.core.coral.models.henryslaw.CoralHenrysLaw;
import insilico.core.coral.models.loael.CoralLOAEL;
import insilico.core.coral.models.loael.general.CoralGeneralLoael;
import insilico.core.coral.models.loael.liver.CoralLiverLoaelNew;
import insilico.core.coral.models.loael.liver.CoralLiverLoel;
import insilico.core.coral.models.noael.general.CoralGeneralNoael;
import insilico.core.coral.models.noael.liver.CoralLiverNoaelNew;
import insilico.core.coral.models.skin.irritation.CoralSkinIrritation;
import insilico.core.coral.models.zebrafish.CoralZebraFish;
import insilico.core.molecule.InsilicoMolecule;
import insilico.core.molecule.conversion.SmilesMolecule;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

public class TestModel {

    public static void main(String[] args) throws Exception {

        ////////////////////

        String sourcefile = "C:\\Users\\alber\\Desktop\\prova\\prova_henry.txt";

        TestCoralModel(sourcefile, new CoralHenrysLaw());
        if (1==1) return;

        ////////////////////

        String[] smi = {
            "CCC([N+]([O-])=O)([N+]([O-])=O)[N+]([O-])=O",
        };

        CoralModel model = new CoralHenrysLaw();

        for (String s : smi) {
            System.out.println("Attributes for molecule: ");
            System.out.println(s);
            System.out.println();
            double v = model.Predict(s);
            System.out.println();
            System.out.println("Predicted value = " + v);
            System.out.println();
            System.out.println("-------------------");
            System.out.println();
        }

    }

    public static void TestCoralModel(String sourcefile, CoralModel model) throws Exception {

        BufferedReader br = new BufferedReader(new FileReader(sourcefile));
        String line;
        int idx = 1;

        while ( (line=br.readLine())!=null) {
            double pred = model.Predict(line);
            System.out.println(idx++ + "\t" + line+"\t" + pred);
        }
    }


}
