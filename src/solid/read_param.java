/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package solid;

import Interface.*;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;

/**
 *
 * @author agung
 */
public class read_param {
    public static void main(String[] args) throws FileNotFoundException, IOException, ClassNotFoundException {
        ObjectInputStream inputStream_ = new ObjectInputStream(new FileInputStream("/home/agung/project/Quantum-20211024T084747Z-001 (2)/jNNDFT_publish_no_gui/JQC_data/ann_new_dft"));

               double krenl[][][][] = (double[][][][]) inputStream_.readObject();
               String bobot="";
               for (int i = 0; i < krenl.length; i++) {
                   //System.out.println("batas "+i);
                    for (int j = 0; j < krenl[i].length; j++) {
                        for (int k = 0; k < krenl[i][j].length; k++) {
                            for (int l = 0; l < krenl[i][j][k].length; l++) {
                                String kk="";
                                if(i==0){
                                kk="w";
                                }
                                else if(i==1){
                                kk="p";
                                } else if(i==2){
                                kk="n";
                                }
                                //System.out.println("$"+kk+"_{"+j+" "+k+" "+l+" }$="+krenl[i][j][k][l]+""); 
                              bobot+=("w("+(i+1)+","+(j+1)+","+(k+1)+","+(l+1)+")="+krenl[i][j][k][l]+""); 
                            }
                        }

                    }

                }
 
    }
   
}
