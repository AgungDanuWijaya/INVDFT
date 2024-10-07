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
import java.util.logging.Level;
import java.util.logging.Logger;
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
        int clus = kernel.a.cluster_num;
        sebar_dist sb = new sebar_dist();

        int mulai = 0;
        System.out.println("run.inversi.run_train()1");

        double best = 9999999;
        double best20 = 9999999;
        double x_new[] = new double[get_le()];
        double del =0.01;//0.075;
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
        dbprocess db3 = new dbprocess(database.local, kernel);
        db3.setData(" DELETE FROM `Quantum`.`error` ");
        db3.closeDB();
        for (int i = 0; i < 20000; i++) {
            int max = 27;
            double save[][] = new double[max][];
            best__ = 9999;

            for (int j = 0; j < max; j++) {
                double in[] = new double[x_new.length];

                for (int k = 0; k < x_new.length; k++) {

                    double x = (x_new[k] - del) + Math.random() * (x_new[k] + del - (x_new[k] - del));
                    if (Math.abs(x_new[k]) < Math.pow(10, -5)) {
                        // x = x_new[k];
                    }
                    if (j == 0 & i == 0) {
                        x = x_new[k];
                    }
                    in[k] = x;
                }

                double w__[][][][] = convert(in);

                save_monte(w__);
                print(convert(in));
                System.out.println("run.inversi_monte.run_train()");
                double tot = 0;
                System.out.println("run.inversi_monte.run_train()sdfsfsdfs");

                int num_cluster = 1;
                dbprocess db2 = new dbprocess(database.local, kernel);
                db2.setData("DELETE FROM `Quantum`.`cluster_monte`");
                System.out.println("run.inversi_monte.run_train()sdfsfsdfs");
                for (int k = 0; k < num_cluster; k++) {
                    String hjhj = "INSERT INTO `Quantum`.`cluster_monte` (`cluster`) VALUES ('" + (k) + "');";
                    db2.setData(hjhj);
                }

                db2.closeDB();

                int y = 1;
                while (y == 1) {
                    dbprocess db1 = new dbprocess(database.local, kernel);
                    String jum = db1.getSingleData("SELECT count(cluster) FROM Quantum.cluster_monte where eror is not null");
                    if (jum.equals(num_cluster + "")) {
                        y = 3;
                        String er = db1.getSingleData("SELECT sum(eror) FROM Quantum.cluster_monte;");
                        tot = Double.parseDouble(er);
                    }

                    db1.closeDB();

                }
                System.out.println(tot + "  sdas");

                save[j] = in;
                if (Math.abs(tot) < best) {
                    best__ = j;
                    best = Math.abs(tot);

                    String Q = "INSERT INTO `Quantum`.`error` (`from1,`errorcol`) VALUES ('monte','" + best + "');";
                    Q = "INSERT INTO `Quantum`.`error` (`errorcol`, `from`) VALUES ('" + best + "', 'monte');";
                    dbprocess db = new dbprocess(database.local, kernel);
                    db.setData(Q);
                    db.closeDB();

                }
                System.err.println(del + "  " + i + " " + j + "================" + tot);

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
/*this.kernel.weight[0][0][ 0][ 0]=-0.3756732069207656;
this.kernel.weight[0][0][ 0][ 1]=0.0372864439641959;
this.kernel.weight[0][0][ 0][ 2]=-0.0535584357815186;
this.kernel.weight[0][0][ 0][ 3]=0.06001283232191324;
this.kernel.weight[0][0][ 1][ 0]=-0.8435838546820698;
this.kernel.weight[0][0][ 1][ 1]=0.016583254533489526;
this.kernel.weight[0][0][ 1][ 2]=-0.575859144650212;
this.kernel.weight[0][0][ 1][ 3]=0.11302410777328208;
this.kernel.weight[0][1][ 0][ 0]=-0.6452521087298986;
this.kernel.weight[0][1][ 0][ 1]=-0.19309300009418554;
this.kernel.weight[0][1][ 0][ 2]=1.4843288394899405;
this.kernel.weight[0][1][ 0][ 3]=0.07574268004520049;
this.kernel.weight[0][1][ 1][ 0]=-0.4114262737886443;
this.kernel.weight[0][1][ 1][ 1]=-1.0086136115825157;
this.kernel.weight[0][1][ 1][ 2]=0.5897572801893132;
this.kernel.weight[0][1][ 1][ 3]=-0.3740198361545266;
this.kernel.weight[0][1][ 2][ 0]=-0.4242404492499216;
this.kernel.weight[0][1][ 2][ 1]=-0.8018209869746168;
this.kernel.weight[0][1][ 2][ 2]=-0.2853684259514026;
this.kernel.weight[0][1][ 2][ 3]=-0.39737803787212256;
this.kernel.weight[0][1][ 3][ 0]=0.4398071224784345;
this.kernel.weight[0][1][ 3][ 1]=0.8063886571523703;
this.kernel.weight[0][1][ 3][ 2]=0.8126767757711494;
this.kernel.weight[0][1][ 3][ 3]=0.1967935436684809;
this.kernel.weight[0][2][ 0][ 0]=-0.2928371529296105;
this.kernel.weight[0][2][ 1][ 0]=0.04519555053862858;
this.kernel.weight[0][2][ 2][ 0]=-0.3007688471795118;
this.kernel.weight[0][2][ 3][ 0]=0.15242612042092238;
this.kernel.weight[1][0][ 0][ 0]=1.3704605815838402;
this.kernel.weight[1][0][ 0][ 1]=1.189435179087346;
this.kernel.weight[1][0][ 0][ 2]=1.4661796137300809;
this.kernel.weight[1][0][ 0][ 3]=1.5782596917046177;
this.kernel.weight[1][0][ 1][ 0]=1.131112451822403;
this.kernel.weight[1][0][ 1][ 1]=0.8779688736668849;
this.kernel.weight[1][0][ 1][ 2]=1.18689469144575;
this.kernel.weight[1][0][ 1][ 3]=0.9787456410014665;
this.kernel.weight[1][1][ 0][ 0]=0.9062977910133232;
this.kernel.weight[1][1][ 0][ 1]=0.7831169888584903;
this.kernel.weight[1][1][ 0][ 2]=1.014322387736647;
this.kernel.weight[1][1][ 0][ 3]=0.6986785550822235;
this.kernel.weight[1][1][ 1][ 0]=1.5578639731005586;
this.kernel.weight[1][1][ 1][ 1]=1.0649231266177421;
this.kernel.weight[1][1][ 1][ 2]=1.9376093285731484;
this.kernel.weight[1][1][ 1][ 3]=0.9513864567362702;
this.kernel.weight[1][1][ 2][ 0]=0.9185972280430644;
this.kernel.weight[1][1][ 2][ 1]=1.297650490285165;
this.kernel.weight[1][1][ 2][ 2]=1.1301444763981567;
this.kernel.weight[1][1][ 2][ 3]=1.2302756315082892;
this.kernel.weight[1][1][ 3][ 0]=0.7333505803363526;
this.kernel.weight[1][1][ 3][ 1]=0.8938656620841475;
this.kernel.weight[1][1][ 3][ 2]=1.1256654266917872;
this.kernel.weight[1][1][ 3][ 3]=1.4070560759436135;
this.kernel.weight[1][2][ 0][ 0]=1.0057688380933552;
this.kernel.weight[1][2][ 1][ 0]=1.0824009934225478;
this.kernel.weight[1][2][ 2][ 0]=0.8403262861900367;
this.kernel.weight[1][2][ 3][ 0]=0.9440992093515678;
this.kernel.weight[2][0][ 0][ 0]=1.0754671314225182;
this.kernel.weight[2][0][ 0][ 1]=-0.05085552270693653;*/
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
            System.err.println("i7u6y567i89" + ex);
        } catch (IOException ex) {
            System.err.println(ex);
        }
    }

    public void save_monte(double w[][][][]) throws IOException {
        System.out.println("run.inversi_monte.save_monte()");
        ObjectOutputStream outputStream;
        try {
            outputStream = new ObjectOutputStream(new FileOutputStream(kernel.URL_ANN + "_monte"));
            outputStream.writeObject(w);
        } catch (FileNotFoundException ex) {
            System.err.println("dgajds" + ex);
        } catch (IOException ex) {
            System.err.println("etrwef" + ex);
        }///home/agung/project/Quantum-20211024T084747Z-001 (2)/jNNDFT_publish_no_gui/JQC_data
        ////home/agung/ase/jNNDFT_publish_no_gui/JQC_data/ann_new_dft
        sebar_dist sb = new sebar_dist();
        try {
            sb.main1();
        } catch (JSchException ex) {
            Logger.getLogger(inversi.class.getName()).log(Level.SEVERE, null, ex);
        } catch (SftpException ex) {
            Logger.getLogger(inversi.class.getName()).log(Level.SEVERE, null, ex);
        }
          System.out.println("run.inversi_monte.save_monte()");
    }

}
