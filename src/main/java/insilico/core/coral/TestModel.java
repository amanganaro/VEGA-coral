package insilico.core.coral;

import insilico.core.coral.models.bloodbrain.CoralBloodBrainBarrier;
import insilico.core.coral.models.noec.harlequin.CoralNoecHarlequin;

import java.io.BufferedReader;
import java.io.FileReader;

public class TestModel {

    public static void main(String[] args) throws Exception {

        ////////////////////

        String sourcefile = "C:\\Users\\alber\\Desktop\\avian.txt";

        TestCoralModel(sourcefile, new CoralBloodBrainBarrier());
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
