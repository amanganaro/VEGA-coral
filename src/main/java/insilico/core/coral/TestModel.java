package insilico.core.coral;

import insilico.core.coral.models.bloodbrain.CoralBloodBrainBarrier;
import insilico.core.coral.models.noael.CoralNoaelRepeatedRat;
import insilico.core.coral.models.noec.avianrepro.CoralNoecAvianRepro;
import insilico.core.coral.models.noec.harlequin.CoralNoecHarlequin;
import insilico.core.coral.models.retentiontime.CoralRetentionPesticides;

import java.io.BufferedReader;
import java.io.FileReader;

public class TestModel {

    public static void main(String[] args) throws Exception {

        ////////////////////

        String sourcefile = "C:\\Users\\alber\\Desktop\\ret.txt";

        TestCoralModel(sourcefile, new CoralRetentionPesticides());
        if (1==1) return;

        ////////////////////

        String[] smi = {
            "O=[N+]([O-])c1cc(cc(c1N(CCC)CCC)[N+](=O)[O-])S(=O)(=O)N",
        };

        CoralModel model = new CoralNoecAvianRepro();

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
