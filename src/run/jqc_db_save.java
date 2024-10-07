package run;

import Interface.Interface;
import defacom.db_cmd;
import function.main_function;
import geo_molecule.data_geo;
import org.json.simple.parser.ParseException;

/**
 *
 * @author Agung Danu Wijaya
 */
public class jqc_db_save {

    public void save(String name_molecule, Interface s) throws ClassNotFoundException, InterruptedException {
        main_function kernel = new main_function();
        kernel.a = s;
        db_cmd j = new db_cmd();
        try {
            j.get_mole(name_molecule, kernel);
        } catch (ParseException ex) {
            System.out.println(ex);
        }
        double[][] R = j.R;
        int numE[] = j.numE;
        String atom[] = new String[numE.length];
        for (int i = 0; i < atom.length; i++) {
            atom[i] = kernel.atom[numE[i] - 1];
        }
        data_geo mole = new data_geo(numE, atom, kernel.mp.mdot(R, kernel.geo.bondLenght));
        kernel.geo.data.put(name_molecule, mole);
        kernel.tipebasis = "6-31g";
        kernel.verbose = 0;
        kernel.URL_int = kernel.a.base + "/int/" + name_molecule + "";
        driver_abinitio drv = new driver_abinitio();
        drv.RHFsetdata(name_molecule, kernel.URL_int, kernel);
    }
}
