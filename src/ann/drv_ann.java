package ann;

import ann.ann;
import com.sun.jna.Native;
import function.main_function;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.logging.Level;
import java.util.logging.Logger;
import jqc.CInterface;

public class drv_ann {

    String libName = "/home/agung/project/Quantum-20211024T084747Z-001 (2)/jNNDFT_publish_no_gui/lda.so";

    public double[] Ex(main_function kernel, double input[]) {
        double Ex_ann[] = new double[input.length];
        ann c = new ann();
        for (int i = 0; i < input.length; i++) {
            double[] rho = {input[i]};
            Ex_ann[i] = kernel.weight[2][0][0][1] * Math.pow(input[i], kernel.weight[2][0][0][0]);
            //Ex_ann[i] = c.f_ANN(rho, kernel.c_node, kernel.weight, kernel.thv, kernel.th);
        }
        return Ex_ann;
    }

    double asinh(double x) {
        return Math.log(x + Math.sqrt(x * x + 1.0));
    }

    public double[] Exc__(main_function kernel, double input[], double gama[]) {
        ObjectInputStream inputStream_ = null;
        double krenl[][][][] = null;
        try {
            inputStream_ = new ObjectInputStream(new FileInputStream("/home/agung/project/Quantum-20211024T084747Z-001 (2)/jNNDFT_publish_no_gui/JQC_data/ann_new_dft_1"));
        } catch (FileNotFoundException ex) {
            Logger.getLogger(drv_ann.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(drv_ann.class.getName()).log(Level.SEVERE, null, ex);
        }
        try {
            krenl = (double[][][][]) inputStream_.readObject();
        } catch (IOException ex) {
            Logger.getLogger(drv_ann.class.getName()).log(Level.SEVERE, null, ex);
        } catch (ClassNotFoundException ex) {
            Logger.getLogger(drv_ann.class.getName()).log(Level.SEVERE, null, ex);
        }

        double Ex_ann[] = new double[input.length];
        ann c = new ann();
        int c_node[] = {2, 4, 4, 1};
        for (int i = 0; i < input.length; i++) {
            double[] rho = {input[i], Math.pow(input[i], krenl[2][0][0][0]) * Math.pow(gama[i], krenl[2][0][0][1])};

            if (input[i] < Math.pow(10, -11)) {
                double[] o = {0};

                rho = o;
            } else if (gama[i] < Math.pow(10, -11)) {
                double[] o = {0};

                rho = o;
            } else if (input[i] < Math.pow(10, -11) & gama[i] < Math.pow(10, -11)) {
                double[] o = {0};

                rho = o;
            }

            double rt = c.f_ANN(rho, kernel.c_node, kernel.weight, kernel.thv, kernel.th);
            if (Double.isNaN(rt)) {
                rt = 1E-10;
            }
            if (Double.isInfinite(rt)) {
                rt = 100;
            }

            Ex_ann[i] = rt;
        }
        return Ex_ann;

    }

    public double[] Exc(main_function kernel, double input[], double gama[]) {
        double Ex_ann[] = new double[input.length];
        ann c = new ann();
        for (int i = 0; i < input.length; i++) {
            double[] rho = {input[i], Math.pow(input[i], kernel.weight[2][0][0][0]) * Math.pow(gama[i], kernel.weight[2][0][0][1])};

            if (input[i] < Math.pow(10, -11)) {
                double[] o = {0, 0};

                rho = o;
            } else if (gama[i] < Math.pow(10, -11)) {
                double[] o = {0, 0};

                rho = o;
            } else if (input[i] < Math.pow(10, -11) & gama[i] < Math.pow(10, -11)) {
                double[] o = {0, 0};

                rho = o;
            }
            double rt = c.f_ANN(rho, kernel.c_node, kernel.weight, kernel.thv, kernel.th);
            if (Double.isNaN(rt)) {
                rt = 1E-10;
            }
            if (Double.isInfinite(rt)) {
                rt = 100;
            }

            Ex_ann[i] = rt;
            //System.out.println(input[i]+" ,"+gama[i]+" ,"+rt);

        }
        //double Ex_1[]=Exc__(kernel, input, gama);
        //Ex_ann=kernel.mp.mdot(Ex_ann,0.3);
        //Ex_1=kernel.mp.mdot(Ex_1, 0.7);
        //Ex_ann=kernel.mp.adddot(Ex_ann, Ex_1);
        return Ex_ann;
    }

    public double[] Exc_cor__(int in, main_function kernel, double input[], double input_b[], double gama[], double gama_b[], double gama_nn[]) {
        double Ex_ann[] = new double[input.length];
        ann c = new ann();
        for (int i = 0; i < input.length; i++) {

            double a = input[i] + 1E-20;
            double b = input_b[i] + 1E-20;
            double gaa_ = gama[i];
            double gbb_ = gama_b[i];

            double n = a + b;
            double gnn_ = gama_nn[i]; // double gnn_ =  gaa_+gbb_;

            double[] rho = {a, b};
            if (a < Math.pow(10, -11) || b < Math.pow(10, -11)) {
                double[] o = {0, 0};

                rho = o;
            }
         double  rt = c.f_ANN(rho, kernel.c_node, kernel.weight, kernel.thv, kernel.th);
            if (Double.isNaN(rt)) {
                rt = 1E-10;
            }
            if (Double.isInfinite(rt)) {
                rt = 100;
            }

            Ex_ann[i] = rt;
        }

        return Ex_ann;
    }

    public double[] Exc_cor(int in, main_function kernel, double input[], double input_b[], double gama[], double gama_b[], double gama_nn[]) {
        double Ex_ann[] = new double[input.length];
        CInterface demo = (CInterface) Native.load(libName, CInterface.class);

        for (int i = 0; i < input.length; i++) {
            double A = 0.04918;
            double B = 0.132;
            double C = 0.2533;
            double Dd = 0.349;

            double CF = (3.0 / 10.0) * Math.pow(3.0 * Math.pow(Math.PI, 2), 2.0 / 3.0);

            double gamma_s_[] = {gama[i], gama_b[i]};

            double rho__[] = {input[i] + 1E-20, input_b[i] + 1E-20};

            double a = rho__[0];
            double b = rho__[1];
            double gaa_ = gamma_s_[0];
            double gbb_ = gamma_s_[1];

            double n = a + b;
            double gnn_ = gama_nn[i]; // double gnn_ =  gaa_+gbb_;

            // double gnn_ =  gaa_+gbb_;
            double gaa = Math.pow(gaa_, 1);
            double gbb = Math.pow(gbb_, 1);
            double gnn = Math.pow(gnn_, 1);

            double icbrtn = Math.pow(n, -1.0 / 3.0);
            double P = 1 / (1 + Dd * icbrtn);
            double omega = Math.exp(-C * icbrtn) * P * Math.pow(n, -11.0 / 3.0);
            double delta = icbrtn * (C + Dd * P);
            double n2 = n * n;
            double f = -A * (4 * a * b * P / n
                    + B * omega
                    * (a * b
                    * (Math.pow(2, 11.0 / 3.0) * CF
                    * (Math.pow(a, 8.0 / 3.0) + Math.pow(b, 8.0 / 3.0))
                    + (47.0 - 7.0 * delta) * gnn / 18.0
                    - (2.5 - delta / 18.0) * (gaa + gbb)
                    - (delta - 11.0) / 9.0 * (a * gaa + b * gbb) / n)
                    - 2.0 / 3.0 * n2 * gnn + (2.0 / 3.0 * n2 - a * a) * gbb
                    + (2.0 / 3.0 * n2 - b * b) * gaa));
            Ex_ann[i] = f;

            //demo.gga_x(rho__, gamma_s_, Fx_, dFxdn_, dero_, Fx_.length, 131);
            //Ex_ann[i]=Fx_[0];
            //System.out.println(a+", "+b+", "+","+gaa_+","+gbb_+"," +Ex_ann[i]+", "+Fx_[0]+" ,"+f);
        }

        return Ex_ann;
    }

    public double[] Exc_cor_(int in, main_function kernel, double input[], double input_b[], double gama[], double gama_b[], double gama_nn[]) {
        double Ex_ann[] = new double[input.length];
        CInterface demo = (CInterface) Native.load(libName, CInterface.class);
        for (int i = 0; i < input.length; i++) {
            /*
             */
            double A = kernel.weight[0][0][0][0];
            double B = kernel.weight[1][0][0][0];
            double C = kernel.weight[2][0][0][0];
            double Dd = kernel.weight[2][0][0][1];
            //System.out.println(Dd);
            
            
            

            double CF = (3.0 / 10.0) * Math.pow(3.0 * Math.pow(Math.PI, 2), 2.0 / 3.0);

            double gamma_s_[] = {gama[i], gama_b[i]};

            double rho__[] = {input[i] + 1E-20, input_b[i] + 1E-20};

            double a = rho__[0];
            double b = rho__[1];
            double gaa_ = gamma_s_[0];
            double gbb_ = gamma_s_[1];

            double n = a + b;
            double gnn_ = gama_nn[i]; // double gnn_ =  gaa_+gbb_;

            // double gnn_ =  gaa_+gbb_;
            double gaa = Math.pow(gaa_, 1);
            double gbb = Math.pow(gbb_, 1);
            double gnn = Math.pow(gnn_, 1);

            double icbrtn = Math.pow(n, -1.0 / 3.0);
            double P = 1 / (1 + Dd * icbrtn);
            double omega = Math.exp(-Math.abs(C) * icbrtn) * P * Math.pow(n, -11.0 / 3.0);
            double delta = icbrtn * (C + Dd * P);
            double n2 = n * n;
            double f = -A * (4 * a * b * P / n
                    + B * omega
                    * (a * b
                    * (Math.pow(2, 11.0 / 3.0) * CF
                    * (Math.pow(a, 8.0 / 3.0) + Math.pow(b, 8.0 / 3.0))
                    + (47.0 - 7.0 * delta) * gnn / 18.0
                    - (2.5 - delta / 18.0) * (gaa + gbb)
                    - (delta - 11.0) / 9.0 * (a * gaa + b * gbb) / n)
                    - 2.0 / 3.0 * n2 * gnn + (2.0 / 3.0 * n2 - a * a) * gbb
                    + (2.0 / 3.0 * n2 - b * b) * gaa));
            Ex_ann[i] = f;
            //System.out.println(a+" "+b+" "+gaa+" "+gbb+" "+gnn+" "+f);

            //demo.gga_x(rho__, gamma_s_, Fx_, dFxdn_, dero_, Fx_.length, 131);
            //Ex_ann[i]=Fx_[0];
            //System.out.println(a+", "+b+", "+","+gaa_+","+gbb_+"," +Ex_ann[i]+", "+Fx_[0]+" ,"+f);
        }

        return Ex_ann;
    }

    public double[] Exc_meta(main_function kernel, double input[], double gama[], double gamas[]) {
        double Ex_ann[] = new double[input.length];
        ann c = new ann();
        for (int i = 0; i < input.length; i++) {
            double[] rho = {Math.pow(input[i], kernel.weight[2][0][0][0]) * Math.pow(gama[i], kernel.weight[2][0][0][1]), Math.pow(input[i], kernel.weight[2][0][0][2]) * Math.pow(gamas[i], kernel.weight[2][0][0][3]), input[i]};

            if (input[i] < Math.pow(10, -11)) {
                double[] o = {0, 0, 0};

                rho = o;
            } else if (gama[i] < Math.pow(10, -11)) {
                double[] o = {0, 0, 0};

                rho = o;
            } else if (input[i] < Math.pow(10, -11) & gama[i] < Math.pow(10, -11)) {
                double[] o = {0, 0, 0};

                rho = o;
            }

            Ex_ann[i] = c.f_ANN(rho, kernel.c_node, kernel.weight, kernel.thv, kernel.th);

        }
        return Ex_ann;
    }
}
