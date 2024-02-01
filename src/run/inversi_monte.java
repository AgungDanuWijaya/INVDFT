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
import jqc.randy;
import static run.monte_carlo.y;

public class inversi_monte {

    int ak = 0;
    int aj = 0;
    public HashMap<Integer, Double> jacobian = new HashMap<Integer, Double>();
    Interface a;
    main_function kernel;

    public inversi_monte(Interface a) {
        this.a = a;
        main_function kernel = new main_function();
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

    public void test_all(int simpul[], String name_ann) throws ClassNotFoundException, IOException {
        prepare_cmn(simpul, name_ann);
        double r = 0;
        data a = new data();
        String names[] = a.ae.split("}");
        for (int i = 0; i < names.length; i++) {
            String name[] = names[i].split(",");
            double en[] = new double[name.length];

            r = hitung(name[0]);
            en[0] = r;
            for (int j = 1; j < name.length; j++) {
                double e = hitung(name[j]);
                en[j] = e;
                r -= e;
            }
            System.out.println(name[0] + " : " + r + " ");
        }
    }

    public void TE(int simpul[], String name_ann) throws ClassNotFoundException, IOException {
        prepare_cmn(simpul, name_ann);
        double r = 0;
        data a = new data();
        String names[] = a.ae.split("}");
        for (int i = 0; i < names.length; i++) {
            String name[] = names[i].split(",");
            double en[] = new double[name.length];
            System.out.print(name[0] + ",");
            r = hitung(name[0]);

        }
    }

    public void prepare_cmn(int simpul[], String name) throws ClassNotFoundException {
        this.kernel.a = this.a;
        this.kernel.tipebasis = a.basis;
        this.kernel.verbose = 1; // 1 tambilkan proses
        String j[] = kernel.a.ann_conf.split(",");
        //int simpul[] = {1, 2, 2, 1};
        this.kernel.c_node = simpul;
        this.kernel.c_node = simpul;
        this.kernel.URL_ANN = kernel.a.base + "/" + name;
        // this.kernel.URL_ANN = kernel.a.base.getText() + "/" + "ann_new_dft";
        //this.kernel.Ex = "A-metaNN";
         this.kernel.Ex = this.kernel.a.exc_tipe;
        if(this.kernel.a.user_param.equals("1")){
        this.kernel.weight = new double[3][][][];
        }else{
         this.kernel.weight = new double[2][][][];
        }
        init(1);
    }

    public double[][][][] convert(double[] par) {
        double w____[][][][] = new double[this.kernel.weight.length][][][];
        int tan = 0;
        for (int i = 0; i < this.kernel.weight.length; i++) {
            double w___[][][] = new double[this.kernel.weight[i].length][][];
            for (int j = 0; j < this.kernel.weight[i].length; j++) {
                double w__[][] = new double[this.kernel.weight[i][j].length][];
                for (int k = 0; k < this.kernel.weight[i][j].length; k++) {
                    double w_[] = new double[this.kernel.weight[i][j][k].length];

                    for (int m = 0; m < this.kernel.weight[i][j][k].length; m++) {
                        w_[m] = par[tan];
                        tan++;

                    }
                    w__[k] = w_;
                }
                w___[j] = w__;
            }
            w____[i] = w___;
        }
        return w____;
    }

    public void print(double w[][][][]) {
        for (int i = 0; i < this.kernel.weight.length; i++) {
            for (int j = 0; j < this.kernel.weight[i].length; j++) {
                for (int k = 0; k < this.kernel.weight[i][j].length; k++) {
                    for (int m = 0; m < this.kernel.weight[i][j][k].length; m++) {

                        try {
                            System.out.print((w[i][j][k][m] + "") + " : ");
                        } catch (Exception e) {
                            System.out.print(w[i][j][k][m] + "" + " : ");
                        }

                    }
                }
            }
        }
    }

    public void run_train(String name[], double re[][], int simpul[], String name_ann) throws ClassNotFoundException, IOException {
        prepare_cmn(simpul, name_ann);
        /*    String name[] = {"H2S,S,H,H}", "Na2,Na,Na}", "H2O2,H,H,O,O}", "HCl,H,Cl}", "HF,H,F}", "Cl2,Cl,Cl}", "H2,H,H}",
            "LiF,Li,F}", "LiH,Li,H}", "CH,C,H}", "OH,O,H}", "Cl2"};
         */

        double re_t[][] = new double[re.length][re[0].length];
        double delta = Math.pow(10, -7);
        double learning = 0.9;//0.9
        double peredam = 0.1;
        double iter = 10000;
        System.out.println("run.inversi.run_train()1");
        init(1);
        save(this.kernel.weight);
        int mulai = 0;
        System.out.println("run.inversi.run_train()1");

        double best = 9999999;
        double x_new[] = new double[28];
        double del = 0.02;//0.016739573143016696
        randy rn = new randy();
        rn.init();

        int tan = 0;
        for (int i = 0; i < this.kernel.weight.length; i++) {
            for (int j = mulai; j < this.kernel.weight[i].length; j++) {
                for (int k = 0; k < this.kernel.weight[i][j].length; k++) {
                    for (int m = 0; m < this.kernel.weight[i][j][k].length; m++) {
                        int ijk[] = {i, j, k, m};
                        try {
                            System.out.print((this.kernel.weight[i][j][k][m] + "") + " : ");
                        } catch (Exception e) {
                            System.out.print(this.kernel.weight[i][j][k][m] + "" + " : ");
                        }

                    }
                }
            }
        }
        for (int i = 0; i < this.kernel.weight.length; i++) {
            for (int j = mulai; j < this.kernel.weight[i].length; j++) {
                for (int k = 0; k < this.kernel.weight[i][j].length; k++) {
                    for (int m = 0; m < this.kernel.weight[i][j][k].length; m++) {
                        x_new[tan] = this.kernel.weight[i][j][k][m];
                        tan++;
                    }
                }
            }
        }

        double we_Re[][][][] = convert(x_new);
        tan = 0;
        int best__ = 0;
        for (int i = 0; i < 20; i++) {
            double save[][] = new double[100][];
            best__ = 9999;
            for (int j = 0; j < 30; j++) {
                double in[] = new double[x_new.length];

                for (int k = 0; k < x_new.length; k++) {

                    double x = (x_new[k] - del) + Math.random() * (x_new[k] + del - (x_new[k] - del));
                    if (Math.abs(x_new[k]) < Math.pow(10, -5)) {
                        x = x_new[k];
                    }
                    if (j == 0 & i == 0) {
                        x = x_new[k];
                    }
                    in[k] = x;
                }

                double w__[][][][] = convert(in);
                this.kernel.weight = w__;
                double y_ = y(in);

                double tot = 0;
                print(convert(in));
                System.out.println("");
                for (int l = 0; l < name.length; l++) {
                    double grad = grad(name[l], w__);

                    System.out.println(name[l] + " " + re[l][0] + "  " + grad);
                    tot += Math.abs((re[l][0] - grad));
                    if (tot > best) {
                        l = 99999;
                    }
                }
                System.err.println(del + "  " + i + " " + j + "================" + tot);
                save[j] = in;
                if (Math.abs(tot) < best) {
                    best__ = j;
                    best = Math.abs(tot);
                    for (int k = 0; k < in.length; k++) {
                      //  x_new[k] = in[k];
                    }

                    String Q = "INSERT INTO `Quantum`.`error` (`from1,`errorcol`) VALUES ('monte','" + best + "');";
                    Q = "INSERT INTO `Quantum`.`error` (`errorcol`, `from`) VALUES ('" + best + "', 'monte');";
                    dbprocess db = new dbprocess(database.local, kernel);
                    db.setData(Q);
                    db.closeDB();

                }
            }
            if (best__ != 9999) {
                x_new = save[best__];
                save(convert(x_new));
            }
            del /= 1.08;
            System.out.println(best);
            print(convert(x_new));

        }
        save(convert(x_new));

    }

    public double[][] a_t_a(double a[][], double b[][]) {
        double r[][] = new double[a[0].length][a.length];
        for (int i = 0; i < a[0].length; i++) {
            for (int j = 0; j < a.length; j++) {
                r[i][j] = a[j][i];
            }
        }
        double dot[][] = new double[r.length][b[0].length];
        for (int i = 0; i < r.length; i++) {
            for (int j = 0; j < b[0].length; j++) {
                for (int j2 = 0; j2 < b.length; j2++) {
                    dot[i][j] += r[i][j2] * b[j2][j];
                }
            }
        }
        return dot;
    }

    public synchronized void input(double grad, int tan) {
        this.jacobian.put(tan, grad);
        aj++;
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

    public double grad(String name, double w[][][][]) throws ClassNotFoundException, IOException {
        double t = 0;
        if (name.contains("}") == false) {
            t = hitung(name);
        } else {
            t = ae(name);
        }

        return t;
    }

    public void init(int status) throws ClassNotFoundException {
        int[] node = new int[kernel.c_node.length - 1];
        for (int i = 1; i < kernel.c_node.length; i++) {
            node[i - 1] = kernel.c_node[i];
        }
        double w[][][] = new double[node.length][][];
        double w_p[][][] = new double[node.length][][];

        for (int i = 0; i < w.length; i++) {
            double wel[][] = new double[kernel.c_node[i]][];
            double wel_p[][] = new double[kernel.c_node[i]][];
            for (int j = 0; j < wel.length; j++) {
                double wdummy[] = new double[node[i]];
                wel[j] = this.kernel.mp.adddot(wdummy, 0.0001);
                wel_p[j] = this.kernel.mp.adddot(wdummy, 1.000); // optimise disini
                //wel[j] = this.kernel.mp.adddot(wdummy, 0.01 * Math.random());
                //wel_p[j] = this.kernel.mp.adddot(wdummy, (0.1 * Math.random() + 0.8)); // optimise disini

            }
            w[i] = wel;
            w_p[i] = wel_p;
        }

       double pa_ne[][][] = new double[1][1][];
        pa_ne[0][0] = this.kernel.a.pa_ne;
        

        this.kernel.weight[0] = w;
        this.kernel.weight[1] = w_p;
        if(this.kernel.a.user_param.equals("1")){
        this.kernel.weight[2] = pa_ne;
        }
        if (status == 1) {
            try {
                ObjectInputStream inputStream = new ObjectInputStream(new FileInputStream(kernel.URL_ANN));

                double krenl[][][][] = (double[][][][]) inputStream.readObject();

                for (int i = 0; i < krenl.length; i++) {
                    for (int j = 0; j < krenl[i].length; j++) {
                        for (int k = 0; k < krenl[i][j].length; k++) {
                            for (int l = 0; l < krenl[i][j][k].length; l++) {
                                this.kernel.weight[i][j][k][l] = krenl[i][j][k][l];
                            }
                        }

                    }

                }

                //this.kernel.weight = (double[][][][]) inputStream.readObject();
                inputStream.close();
            } catch (FileNotFoundException ex) {
            } catch (IOException ex) {
            }
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

    public void save(double w[][][][]) {
        ObjectOutputStream outputStream;
        try {
            outputStream = new ObjectOutputStream(new FileOutputStream(kernel.URL_ANN));
            outputStream.writeObject(w);
        } catch (FileNotFoundException ex) {
            System.err.println(ex);
        } catch (IOException ex) {
            System.err.println(ex);
        }
    }

}
