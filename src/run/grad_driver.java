package run;

import Interface.Interface;
import com.google.gson.Gson;
import java.io.IOException;

/**
 *
 * @author Agung Danu Wijaya
 */
public class grad_driver {

    public void manis(Interface gui) throws ClassNotFoundException, IOException {
        inversi a = new inversi(gui);
        Gson gson = new Gson();
        param_inv init = gson.fromJson(gui.data_geo, param_inv.class);

       a.run_grad(init.name, init.re,init.simpul,init.name_ann);
       //a.read_Save(init.simpul,init.name_ann);
       //a.test_all(init.simpul,init.name_ann);
        //a.TE();
    }

}
