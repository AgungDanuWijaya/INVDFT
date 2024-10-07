package solid;

import java.io.File;  // Import the File class
import java.io.FileNotFoundException;  // Import this class to handle errors
import java.util.Scanner; // Import the Scanner class to read text files

public class baca_band {

    public static double main(String[] args) {
      double en=0;
        try {
            File myObj = new File("/home/agung/project/solid/q-e-qe-6.6/PW/aa");
            Scanner myReader = new Scanner(myObj);
            String hh = "";
            while (myReader.hasNextLine()) {
                String data = myReader.nextLine();
                hh += (data) + "\n";
            }
            hh = (hh.substring(hh.indexOf("##") + 2).split("\n")[0]);
            hh = hh.replace(" ", "#") + "#";
            System.out.println(hh);
            System.err.println(hh.indexOf("879886"));

            String jk[] = hh.split("#");
            double rt[] = {0, 0};
            int in = 0;
            for (int i = 0; i < jk.length; i++) {
                try {
                    rt[in] = Double.parseDouble(jk[i]);
                    in++;
                } catch (Exception e) {
                }

            }

            System.out.println(rt[0] - rt[1]);
en=rt[0] - rt[1];
            myReader.close();

        } catch (FileNotFoundException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
 return en;
    }
   
}
