/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package run;

import com.google.gson.Gson;
import defacom.database;
import defacom.dbprocess;
import function.main_function;
import java.awt.image.Kernel;
import java.util.ArrayList;
import javax.swing.JOptionPane;

/**
 *
 * @author agung
 */
public class tes {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        String dat="H2S,S,H,H}Na2,Na,Na}Si2,Si,Si}P2,P,P}S2,S,S}NaCl,Na,Cl}BeH,H,Be}HCl,H,Cl}HF,H,F}Cl2,Cl,Cl}H2,H,H}LiF,Li,F}LiH,Li,H}CH,C,H}OH,O,H}H2O,O,H,H}O2,O,O}NH,N,H}Li2,Li,Li}CO,C,O}F2,F,F}HCO,H,C,O}H2O2,H,H,O,O}H2CO,H,H,C,O}CH3,C,H,H,H}CH4,C,H,H,H,H}N2,N,N}C2H2,C,C,H,H}CH3OH,O,C,H,H,H,H}NH3,N,H,H,H}HCN,H,C,N}CN,C,N}NH2,N,H,H}C2H6,C,C,H,H,H,H,H,H}N2H4,N,N,H,H,H,H}SiH3,Si,H,H,H}SiH4,Si,H,H,H,H}PH2,P,H,H}PH3,P,H,H,H}CH3SH,C,S,H,H,H,H}SO2,S,O,O}FCl,F,Cl}CH3Cl,C,H,H,H,Cl}AlCl3,Al,Cl,Cl,Cl}";
        String kl[]=dat.split("}");
        for (int i = 0; i < kl.length; i++) {
            System.out.print("\""+kl[i]+"}\",");
        }
         
           
        
    }

}
