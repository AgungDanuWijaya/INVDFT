
package run;

import defacom.db_cmd;
import function.main_function;
import geo_molecule.data_geo;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import org.json.simple.parser.ParseException;

/**
 *
 * @author Agung Danu Wijaya
 */
public class run {

    main_function kernel;

    public run(main_function kernel) {
        this.kernel = kernel;
    }

    public double run(String name_molecule, String code) throws ClassNotFoundException, IOException {
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
        driver_abinitio drv = new driver_abinitio();
        if (code.equals("DFT_ann")) {
            drv.DFT(name_molecule, kernel);
        } else if (code.equals("HF")) {
            drv.HF(name_molecule, kernel);
        }
        return drv.input;
    }

    public double run_scf(String name_molecule, String code) throws ClassNotFoundException, IOException {
        kernel.URL_int = kernel.a.base + "/int/" + name_molecule + "";
        double r = run(name_molecule, code);
        return r;
    }

    public void init() throws ClassNotFoundException {
        int[] node = new int[kernel.c_node.length - 1];
        for (int i = 1; i < kernel.c_node.length; i++) {
            node[i - 1] = kernel.c_node[i];
        }
        try {
            ObjectInputStream inputStream = new ObjectInputStream(new FileInputStream(kernel.URL_ANN));
            kernel.weight = (double[][][][]) inputStream.readObject();
            inputStream.close();
        } catch (FileNotFoundException ex) {
        } catch (IOException ex) {
        }
        double thv[][] = new double[node.length][];
        double th[] = new double[node.length];
        for (int i = 0; i < thv.length; i++) {
            double w_l[] = new double[node[i]];
            w_l = kernel.mp.adddot(w_l, -0.00);
            thv[i] = w_l;
            th[i] = 0.00;
        }
        kernel.th = th;
        kernel.thv = thv;
    }

}
