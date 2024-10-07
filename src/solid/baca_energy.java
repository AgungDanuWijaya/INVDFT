package solid;

import java.io.File;  // Import the File class
import java.io.FileNotFoundException;  // Import this class to handle errors
import java.util.Scanner; // Import the Scanner class to read text files

public class baca_energy {

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
            //hh = (hh.substring(hh.indexOf("total energy"))).substring(0,hh.indexOf("Ry"));
            //hh = hh.replace(" ", "") ;
try {
            hh=hh.substring(hh.indexOf("@@"));
           hh=hh.substring(hh.indexOf("=")+1,hh.indexOf("Ry") ).replace(" ", "");     
            } catch (Exception e) {
            }
           
           
            try {
                en=Double.parseDouble(hh);
            } catch (Exception e) {
            }
            //System.out.println(hj);
            myReader.close();

        } catch (FileNotFoundException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
        return en;

    }
   
}
