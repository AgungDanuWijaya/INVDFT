/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package Interface;

import java.util.logging.Level;
import java.util.logging.Logger;
import run.ann_dft_hf;
import run.jqc_db_save;

/**
 *
 * @author agung
 */
public class save_int {

    public static void main(String[] args) {
        Interface in = new Interface();
        in.int_stat=true;
        in.data_src="MySQL";
        jqc_db_save a = new jqc_db_save();
        String name = in.data_geo.split(";")[1];
        try {
            a.save(name, in);
        } catch (ClassNotFoundException ex) {
            Logger.getLogger(Interface.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(Interface.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
