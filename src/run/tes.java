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
        Gson gson = new Gson();
        param_inv init = new param_inv();
        String json = gson.toJson(init);
        System.out.println(json);
        
         
           
        
    }

}
