package insilico.core.coral;

import insilico.core.coral.models.ec50.CoralEC50;
import insilico.core.coral.models.henryslaw.CoralHenrysLaw;
import insilico.core.coral.models.loael.CoralLOAEL;
import insilico.core.coral.models.loael.general.CoralGeneralLoael;
import insilico.core.coral.models.loael.liver.CoralLiverLoaelNew;
import insilico.core.coral.models.loael.liver.CoralLiverLoel;
import insilico.core.coral.models.noael.general.CoralGeneralNoael;
import insilico.core.coral.models.noael.liver.CoralLiverNoaelNew;
import insilico.core.coral.models.noec.avianrepro.CoralNoecAvianReproduction;
import insilico.core.coral.models.noec.harlequin.CoralNoecHarlequin;
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

        String sourcefile = "C:\\Users\\alber\\Desktop\\avian.txt";

        TestCoralModel(sourcefile, new CoralNoecAvianReproduction());
        if (1==1) return;

        ////////////////////

        String[] smi = {
            "Clc2ccc(C(OCC=C)Cn1ccnc1)c(Cl)c2",
        };

        CoralModel model = new CoralNoecHarlequin();

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
