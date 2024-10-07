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

public class inversi_th {

    int ak = 0;
    int aj = 0;
    public HashMap<Integer, Double> jacobian = new HashMap<Integer, Double>();
    Interface a;
    main_function kernel;

    public inversi_th(Interface a) {
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
        //this.kernel.Ex = "LDA";
        this.kernel.Ex = "ann_metaGGA";
        //this.kernel.weight = new double[2][][][];
        this.kernel.weight = new double[3][][][];
        init(1);
    }

    public void run_train(String name[], double re[][], int simpul[], String name_ann) throws ClassNotFoundException, IOException {
        prepare_cmn(simpul, name_ann);
        /*    String name[] = {"H2S,S,H,H}", "Na2,Na,Na}", "H2O2,H,H,O,O}", "HCl,H,Cl}", "HF,H,F}", "Cl2,Cl,Cl}", "H2,H,H}",
            "LiF,Li,F}", "LiH,Li,H}", "CH,C,H}", "OH,O,H}", "Cl2"};
         */

        double re_t[][] = new double[re.length][re[0].length];
        double delta = Math.pow(10, -3);
        double learning = 0.7;//0.9
        double peredam = 0.2;
        double iter = 10000;
        System.out.println("run.inversi.run_train()1");
        init(1);
        save(this.kernel.weight);
        int mulai = 0;
        System.out.println("run.inversi.run_train()1");
        for (int it = 0; it < iter; it++) {

            System.out.println("run.inversi.run_train()1");
            int tot = 0;
            for (int i = mulai; i < this.kernel.weight.length; i++) {
                for (int j = 0; j < this.kernel.weight[i].length; j++) {
                    for (int k = 0; k < this.kernel.weight[i][j].length; k++) {
                        for (int m = 0; m < this.kernel.weight[i][j][k].length; m++) {
                            int ijk[] = {i, j, k, m};
                            try {
                                System.out.print((this.kernel.weight[i][j][k][m] + "") + " : ");
                            } catch (Exception e) {
                                System.out.print(this.kernel.weight[i][j][k][m] + "" + " : ");
                            }
                            tot++;
                        }
                    }
                }
            }
            System.out.println("run.inversi.run_train()");
            double jacobi[][] = new double[name.length][];

            for (int l = 0; l < name.length; l++) {
                int tan = 0;
                grad grs[] = new grad[tot];
                for (int i = mulai; i < this.kernel.weight.length; i++) {
                    for (int j = 0; j < this.kernel.weight[i].length; j++) {
                        for (int k = 0; k < this.kernel.weight[i][j].length; k++) {
                            for (int m = 0; m < this.kernel.weight[i][j][k].length; m++) {
                                int ijk[] = {i, j, k, m};
                                try {
                                    System.out.print((this.kernel.weight[i][j][k][m] + "").substring(0, 4) + " : ");
                                } catch (Exception e) {
                                    System.out.print(this.kernel.weight[i][j][k][m] + "" + " : ");
                                }
                                grs[tan] = new grad(a, kernel, delta, ijk, name[l], this.kernel.weight);
                                grs[tan].start();
                                tan++;
                            }
                            //System.out.println("");
                        }
                    }
                }

                for (int i = 0; i < grs.length; i++) {
                    int stat = grs[i].se;
                    while (stat == 0) {

                    }
                    double grad = grs[i].grad;
                    grad = grad + grad * 0.001 * Math.random();
                    input(grad, i);
                    System.out.print(i+ " "+grad+" :");
                }

                System.out.println("Selesai " + name[l]);

                double jaco[] = new double[jacobian.size()];
                for (int i = 0; i < jacobian.size(); i++) {
                    jaco[i] = jacobian.get(i);
                }
                jacobi[l] = jaco;
                if (name[l].contains("}") == false) {
                    re_t[l][0] = hitung(name[l]);
                    System.out.println(" " + re_t[l][0] + " " + re[l][0]);
                } else {
                    re_t[l][0] = ae(name[l]);
                    System.out.println(" " + re_t[l][0] + " " + re[l][0]);
                }
            }
            double A[][] = a_t_a(jacobi, jacobi);
            for (int i = 0; i < A.length; i++) {
                A[i][i] += peredam;
            }

            double B[][] = a_t_a(jacobi, kernel.mp.adddot(re, kernel.mp.mdot(re_t, -1)));
            double er = (kernel.mp.sum(kernel.mp.adddotabs(re, kernel.mp.mdot(re_t, -1))));

            System.out.println(er);
            String Q = "INSERT INTO `Quantum`.`error` (`errorcol`) VALUES ('" + er + "');";
            dbprocess db = new dbprocess(database.local, kernel);
            db.setData(Q);
            if (er < 0.27) {
                break;
            }
            db.closeDB();
            Matrix Aj = new Matrix(A);
            Aj = Aj.inverse();
            Aj = Aj.times(new Matrix(B));
            double dp[][] = Aj.getArray();
            int tan = 0;
            for (int i = mulai; i < this.kernel.weight.length; i++) {
                for (int j = 0; j < this.kernel.weight[i].length; j++) {
                    for (int k = 0; k < this.kernel.weight[i][j].length; k++) {
                        for (int l = 0; l < this.kernel.weight[i][j][k].length; l++) {
                            this.kernel.weight[i][j][k][l] += dp[tan][0] * learning;

                            tan += 1;
                        }
                    }
                }
            }
            save(this.kernel.weight);
        }
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

    public double grad(double delta, int[] i, String name, double w[][][][]) throws ClassNotFoundException, IOException {
        double r = w[i[0]][i[1]][i[2]][i[3]];
        double t = 0;
        if (name.contains("}") == false) {
            t = hitung(name);
        } else {
            t = ae(name);
        }
        w[i[0]][i[1]][i[2]][i[3]] += delta;
        double t1 = 0;
        if (name.contains("}") == false) {
            t1 = hitung(name);
        } else {
            t1 = ae(name);
        }
        w[i[0]][i[1]][i[2]][i[3]] = r;
        double grad = (t1 - t) / delta;
        return grad;
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
                wel[j] = this.kernel.mp.adddot(wdummy, 0.001);
                wel_p[j] = this.kernel.mp.adddot(wdummy, 1.000); // optimise disini
                //wel[j] = this.kernel.mp.adddot(wdummy, 0.01 * Math.random());
                //wel_p[j] = this.kernel.mp.adddot(wdummy, (0.1 * Math.random() + 0.8)); // optimise disini

            }
            w[i] = wel;
            w_p[i] = wel_p;
        }

        double pa_ne[][][] = new double[1][1][4];
        pa_ne[0][0][0] = 1.0;
        pa_ne[0][0][1] = 0.000000001;
        pa_ne[0][0][2] = 1.0;
        pa_ne[0][0][3] = 0.0001;

        this.kernel.weight[0] = w;
        this.kernel.weight[1] = w_p;
        //       System.err.println( this.kernel.weight.length+"  jhvhv");
        this.kernel.weight[2] = pa_ne;

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
                ObjectInputStream inputStream_ = new ObjectInputStream(new FileInputStream("/home/agung/project/Quantum-20211024T084747Z-001 (2)/Quantum/JQC_data/simpan/ann_new_dft_GGA"));

                this.kernel.weight_1 = (double[][][][]) inputStream_.readObject();

                //this.kernel.weight = (double[][][][]) inputStream.readObject();
                inputStream_.close();
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

    public void read_Save(int simpul[], String name_ann) throws ClassNotFoundException {
        prepare_cmn(simpul, name_ann);
        String jk = "-0.7126137717842845 : -0.7321365139200837 : 0.014305651529739928 : 0.051349400893782614 : -0.45213408259556237 : -0.4974019544257401 : -0.46367827244109194 : -0.46109027287352194 : -0.7614146985040952 : -0.7723627521806917 : 0.9893799518741702 : 1.290743943246957 : 0.8535045095539872 : 1.2433494307784123 : 1.0083822329706655 : 1.0134642634548847 : 0.9235088144747677 : 0.9776729457158226 : 0.9491227874695594 : 0.9813993999668912 : 1.2289344507265343 : -0.00924223032246411 : 1.1776374433584784 : -0.1622912240218507";
        String hj[] = jk.split(":");
        double krenl[][][][] = this.kernel.weight;
        int in = 0;
        for (int i = 0; i < krenl.length; i++) {
            for (int j = 0; j < krenl[i].length; j++) {
                for (int k = 0; k < krenl[i][j].length; k++) {
                    for (int l = 0; l < krenl[i][j][k].length; l++) {
                        this.kernel.weight[i][j][k][l] = Double.parseDouble(hj[in].replace(" ", ""));
                        in++;
                    }
                }

            }

        }
        save(this.kernel.weight);

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
