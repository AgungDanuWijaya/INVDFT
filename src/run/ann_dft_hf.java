package run;

import Interface.Interface;
import java.util.HashMap;
import function.main_function;
import java.io.IOException;

public class ann_dft_hf {

    Interface a;

    public ann_dft_hf(Interface a) {
        this.a = a;
        main_function kernel = new main_function();
        this.kernel = kernel;
    }
    run main;
    main_function kernel;
    public HashMap<Integer, Double> jacobian = new HashMap<Integer, Double>();

    public void energi() throws ClassNotFoundException, IOException {
        kernel.verbose=a.print;
        kernel.status_int = 1;
        if (this.a.int_stat == false) {
            kernel.status_int = 0;
        }
        kernel.a = this.a;
        run main = new run(kernel);
        this.main = main;
        double r = 0;
        String name = a.data_geo.split(";")[1];
        prepare_cmn();
        kernel.c_mix = Double.parseDouble(kernel.a.mix);
        System.out.println("hh"+kernel.opti);
        r = main.run_scf(name, "DFT_ann_ann");
        System.out.println(r);
    }

    public double energi_dft() throws ClassNotFoundException, IOException {
        kernel.status_int = 1;
          kernel.verbose=a.print;
        if (this.a.int_stat == false) {
            kernel.status_int = 0;
        }
        kernel.a = this.a;
        run main = new run(kernel);
        this.main = main;
        double r = 0;
        String name = a.data_geo.split(";")[1];
        prepare_cmn();
        kernel.c_mix = Double.parseDouble(kernel.a.mix);
        r = main.run_scf(name, "DFT_ann");
       return r;
    }

    public void energi_hf() throws ClassNotFoundException, IOException {
        kernel.status_int = 1;
          kernel.verbose=a.print;
        if (this.a.int_stat== false) {
            kernel.status_int = 0;
        }
        kernel.a = this.a;
        run main = new run(kernel);
        this.main = main;
        double r = 0;
        String name = a.data_geo.split(";")[1];
        prepare_cmn();       
        r = main.run_scf(name, "HF");
        System.out.println(r);
    }

    public void prepare_cmn() throws ClassNotFoundException {
        this.kernel.tipebasis = a.basis;
          kernel.verbose=a.print;
        this.kernel.verbose = 1; // 1 tambilkan proses
        String j[] = kernel.a.ann_conf.split(",");
        int simpul[] = new int[j.length];
        int jk = 0;
        for (Object object : j) {
            simpul[jk++] = Integer.parseInt(object.toString());
        }
        this.kernel.c_node = simpul;
        this.kernel.c_node = simpul;
        this.kernel.URL_ANN = kernel.a.base + "/"+kernel.a.name_exc;
        this.kernel.Ex = a.exc_tipe;
        this.kernel.c_mix = 0.9;
        this.kernel.weight = new double[2][][][];
        main.init();
    }

}
