package run;

import Interface.Interface;
import Jama.Matrix;
import defacom.database;
import defacom.dbprocess;
import java.util.HashMap;
import function.main_function;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.logging.Level;
import java.util.logging.Logger;

public class grad extends Thread {

    Interface a;
    main_function kernel;
    double delta;
    int[] i;
    String name;
    double w[][][][];
    double grad = 0;
    int se = 0;

    public grad(Interface a, main_function kernel, double delta, int[] i, String name, double w[][][][]) {
        this.a = a;
        this.delta = delta;
        this.i = i;
        this.name = name;
        this.w = w;
        this.kernel = kernel;
    }

    public double energi_dft(String name) throws ClassNotFoundException, IOException {
        kernel.status_int = 1;
        if (this.a.int_stat == false) {
            kernel.status_int = 0;
        }
        run main = new run(kernel);
        double r = 0;
        kernel.c_mix = Double.parseDouble(kernel.a.mix);
        r = main.run_scf(name, "DFT_ann");
        //r = main.run_scf(name, "HF");
        return r;
    }

    public double hitung(String name_molecule) throws ClassNotFoundException, IOException {
        kernel.URL_int = "JQC_data/int/" + name_molecule + "";
        double r = energi_dft(name_molecule);
        return r;
    }

    public double ae(String data) throws ClassNotFoundException, IOException {
        double r = 0;
        String name[] = data.replace("}", "").split(",");
        for (int i = 1; i <= 1; i++) {
            double en[] = new double[name.length];
            r = hitung(name[0]);
            en[0] = r;
            for (int j = 1; j < name.length; j++) {
                double e = hitung(name[j]);
                en[j] = e;
                r -= e;
            }
        }
        return r;
    }

    public void run() {
        double r = w[i[0]][i[1]][i[2]][i[3]];
        double t = 0;
        if (name.contains("}") == false) {
            try {
                t = hitung(name);
            } catch (ClassNotFoundException ex) {
                Logger.getLogger(grad.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                Logger.getLogger(grad.class.getName()).log(Level.SEVERE, null, ex);
            }
        } else {
            try {
                t = ae(name);
            } catch (ClassNotFoundException ex) {
                Logger.getLogger(grad.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                Logger.getLogger(grad.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        w[i[0]][i[1]][i[2]][i[3]] += delta;
        double t1 = 0;
        if (name.contains("}") == false) {
            try {
                t1 = hitung(name);
            } catch (ClassNotFoundException ex) {
                Logger.getLogger(grad.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                Logger.getLogger(grad.class.getName()).log(Level.SEVERE, null, ex);
            }
        } else {
            try {
                t1 = ae(name);
            } catch (ClassNotFoundException ex) {
                Logger.getLogger(grad.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                Logger.getLogger(grad.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        w[i[0]][i[1]][i[2]][i[3]] = r;
        double grad = (t1 - t) / delta;
        this.se = 1;
        this.grad = grad;
    }

}
