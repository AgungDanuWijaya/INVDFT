/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package Interface;

import java.util.logging.Level;
import java.util.logging.Logger;
import run.jqc_db_save;

/**
 *
 * @author agung
 */
public class hh {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        String hjj="";
        String jk="\"Cl\",\"He\",\"B,B_}\",\"Ne,Ne_}\",\"Mg,Mg_}\",\"S,S_}\",\"H,H_}\",\"P,P_}\",\"N,N_}\",\"Na,Na_}\",\"Si,Si_}\",\"He,He_}\",\"Be,Be_}\",\"Al,Al_}\",\"Li,Li_}\",\"F,F_}\",\"Cl,Cl_}\",\"C,C_}\",\"O,O_}\",\"H\",\"OH,O,H}\",\"H2,H,H}\",\"H2S,S,H,H}\",\"Na2,Na,Na}\",\"Si2,Si,Si}\",\"P2,P,P}\",\"S2,S,S}\",\"NaCl,Na,Cl}\",\"BeH,H,Be}\",\"HCl,H,Cl}\",\"HF,H,F}\",\"Cl2,Cl,Cl}\",\"H2,H,H}\",\"LiF,Li,F}\",\"LiH,Li,H}\",\"CH,C,H}\",\"OH,O,H}\",\"H2O,O,H,H}\",\"O2,O,O}\",\"NH,N,H}\",\"Li2,Li,Li}\",\"CO,C,O}\",\"F2,F,F}\",\"HCO,H,C,O}\",\"H2O2,H,H,O,O}\",\"H2CO,H,H,C,O}\",\"CH3,C,H,H,H}\",\"CH4,C,H,H,H,H}\",\"N2,N,N}\",\"C2H2,C,C,H,H}\",\"CH3OH,O,C,H,H,H,H}\",\"NH3,N,H,H,H}\",\"HCN,H,C,N}\",\"CN,C,N}\",\"NH2,N,H,H}\",\"C2H6,C,C,H,H,H,H,H,H}\",\"N2H4,N,N,H,H,H,H}\",\"SiH3,Si,H,H,H}\",\"SiH4,Si,H,H,H,H}\",\"PH2,P,H,H}\",\"PH3,P,H,H,H}\",\"CH3SH,C,S,H,H,H,H}\",\"SO2,S,O,O}\",\"FCl,F,Cl}\",\"CH3Cl,C,H,H,H,Cl}\",\"AlCl3,Al,Cl,Cl,Cl}\",\"CH3OCH3,C,H,H,H,O,C,H,H,H}\",\"C3H4,C,C,C,H,H,H,H}\",\"C4H6,C,C,C,C,H,H,H,H,H,H}\",";
    String j[]=jk.replace("\"", "").split("},");
        for (int i = 0; i < j.length; i++) {
            String string[] = j[i].split(",");
            for (int k = 0; k < string.length; k++) {
                String string1 = string[k];
                Interface in = new Interface();
        in.int_stat=true;
        in.data_src="MySQL";
        jqc_db_save a = new jqc_db_save();
        
        in.data_geo=string1+";"+string1+";;";
if(hjj.contains(in.data_geo)==false){        
        System.out.println(in.data_geo);
        String name = in.data_geo.split(";")[1];
        try {
            a.save(name, in);
        } catch (ClassNotFoundException ex) {
            Logger.getLogger(Interface.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(Interface.class.getName()).log(Level.SEVERE, null, ex);
        }
             //   System.out.print(string1+"; "+string1+";;");
 hjj+=in.data_geo;           
}
            System.out.println("");
        }
       
        }
    }
    
}
