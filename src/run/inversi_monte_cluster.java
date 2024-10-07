package run;

import Interface.Interface;
import Interface.sebar_dist;
import Jama.Matrix;
import com.jcraft.jsch.JSchException;
import com.jcraft.jsch.SftpException;
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
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import jqc.randy;
import static run.monte_carlo.y;

public class inversi_monte_cluster {

    int ak = 0;
    int aj = 0;
    public HashMap<Integer, Double> jacobian = new HashMap<Integer, Double>();
    Interface a;
    main_function kernel;

    public inversi_monte_cluster(Interface a) {
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
        if (this.kernel.a.user_param.equals("1")) {
            this.kernel.weight = new double[3][][][];
        } else {
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

    public int get_le() {
        int n = 0;
        for (int i = 0; i < this.kernel.weight.length; i++) {
            for (int j = 0; j < this.kernel.weight[i].length; j++) {
                for (int k = 0; k < this.kernel.weight[i][j].length; k++) {
                    for (int m = 0; m < this.kernel.weight[i][j][k].length; m++) {

                        n++;
                    }
                }
            }
        }
        return n;
    }

    public void run_train(String name[], double re[][], int simpul[], String name_ann, int clus) throws ClassNotFoundException, IOException, InterruptedException {
        prepare_cmn(simpul, name_ann);

        while (true) {
            dbprocess db21 = new dbprocess(database.local, kernel);
            String hjj[] = db21.getRowSetData("SELECT * FROM Quantum.cluster_monte where cluster='" + clus + "' and eror is null;");

            db21.closeDB();
            TimeUnit.SECONDS.sleep(2);
            System.out.println("cekikng data");
            if (hjj != null) {
                init_monte(1);
                print(this.kernel.weight);
                double tot = 0;
                

                for (int l_ = 0; l_ <name.length; l_++) {
                    
                    int l=l_;

                    double grad = grad(name[l], this.kernel.weight);
                    System.out.println();

                    System.out.println(name[l] + " " + re[l][0] + "  " + grad + " " + tot + " " + l);

                    if (l < 2) {
                        tot += Math.abs((re[l][0] - grad) / 30);
                    } else {
                        tot += Math.abs((re[l][0] - grad));
                    }
                    String jj = "UPDATE `Quantum`.`cluster_monte` SET `eror1` = '" + tot + "' WHERE (`cluster` = '" + clus + "');";
                    dbprocess db2 = new dbprocess(database.local, kernel);
                    db2.setData(jj);
                    String bests = db2.getSingleData("SELECT errorcol FROM Quantum.error order by id desc limit 1;");
                    String er = db2.getSingleData("SELECT sum(eror1) FROM Quantum.cluster_monte;");
                    double r = 0;
                    double r1 = 0;
                    try {
                        r = Double.parseDouble(bests);
                    } catch (Exception e) {
                        r=999999;
                    }
                    try {
                        r1 = Double.parseDouble(er);
                    } catch (Exception e) {
                       
                    }
                    if(r1>r){
                    tot=999999;
                        break;
                    }
                    
                   

                    db2.closeDB();
                }
                String jj = "UPDATE `Quantum`.`cluster_monte` SET `eror` = '" + tot + "' WHERE (`cluster` = '" + clus + "');";
                dbprocess db2 = new dbprocess(database.local, kernel);
                db2.setData(jj);
                db2.closeDB();
            }
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

    public double grad(String name, double w[][][][]) throws ClassNotFoundException, IOException {
        double t = 0;
        if (name.contains("}") == false) {
            t = hitung(name);
            if (Double.isNaN(t)) {
                t = 10E9;
            }
        } else {
            t = ae(name);
            if (Double.isNaN(t)) {
                t = 10E9;
            }
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
                wel[j] = this.kernel.mp.adddot(wdummy, this.kernel.a.we);
                wel_p[j] = this.kernel.mp.adddot(wdummy, this.kernel.a.we_p); // optimise disini
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
        if (this.kernel.a.user_param.equals("1")) {
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
                //inputStream.close();
                //ObjectInputStream inputStream_ = new ObjectInputStream(new FileInputStream(this.kernel.a.base+"/"+this.kernel.a.name_exc));
                //this.kernel.weight_1 = (double[][][][]) inputStream_.readObject();
                //this.kernel.weight = (double[][][][]) inputStream.readObject();
                //inputStream_.close();
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

    public void init_monte(int status) throws ClassNotFoundException {
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
                wel[j] = this.kernel.mp.adddot(wdummy, this.kernel.a.we);
                wel_p[j] = this.kernel.mp.adddot(wdummy, this.kernel.a.we_p); // optimise disini
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
        if (this.kernel.a.user_param.equals("1")) {
            this.kernel.weight[2] = pa_ne;
        }

        if (status == 1) {
            try {
                ObjectInputStream inputStream = new ObjectInputStream(new FileInputStream(kernel.URL_ANN + "_monte"));

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
                //inputStream.close();
                //ObjectInputStream inputStream_ = new ObjectInputStream(new FileInputStream(this.kernel.a.base+"/"+this.kernel.a.name_exc));
                //this.kernel.weight_1 = (double[][][][]) inputStream_.readObject();
                //this.kernel.weight = (double[][][][]) inputStream.readObject();
                //inputStream_.close();
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
