/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package Interface;

import com.google.gson.Gson;
import java.io.IOException;
import run.inversi_driver;
import run.inversi_monte;
import run.inversi_monte_cluster;
import run.param_inv;

/**
 *
 * @author agung
 */
public class inversi_cluster {
    public static void main(String[] args) throws ClassNotFoundException, IOException, InterruptedException {
      
        Interface gui=new Interface();
        inversi_monte_cluster a = new inversi_monte_cluster(gui);
        Gson gson = new Gson();
        param_inv init = gson.fromJson(gui.data_geo, param_inv.class);

        a.run_train(init.name, init.re,init.simpul,init.name_ann,1);
    }
}
