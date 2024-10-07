/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jqc;

import com.sun.jna.Native;
import tools.matrix_operation;

/**
 *
 * @author agung
 */
public class NewMain1 {

    /**
     * @param args the command line arguments
     */
    static double asinh(double a) {

        return Math.log(Math.sqrt(a * a + 1.0d) + a);

    }

    public static void main(String[] args) {
        matrix_operation ao = new matrix_operation();
        double rho[] = {0.2};
        double dero[] = new double[rho.length];
        double gamma_s[] = {0.4};

        double rho13[] = ao.powdot(rho, 4.0 / 3.0);

        double x[] = ao.divdot(ao.powdot(gamma_s, 0.5), rho13);
        //  ao.disp(x);
        //x = ao.divdot(x, rho);
        // x[0]=1.0;

        double b = 0.0042;
        double b88_g[] = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            double b88_g_ = -1.5 * Math.pow(3.0 / 4.0 / Math.PI, 1.0 / 3.0) - b * x[i] * x[i] / (1.0 + 6.0 * b * x[i] * asinh(x[i]));
            //    double b88_g_ = -b * x[i] * x[i] / (1.0 + 6.0 * b * x[i] * asinh( Math.pow(1.0, 1.0 / 3.0)*x[i]));
            System.out.println(b88_g_);
            b88_g[i] = b88_g_;
        }

        double Fx[] = ao.mdot(ao.powdot(rho, 4.0 / 3.0), b88_g);
        //ao.disp(Fx);
        // ao.disp(ao.divdot(Fx,rho));
        double[] dFxdn = ao.mdot(rho13, (4.0 / 3.0));

        String libName = "/home/agung/project/Quantum-20211024T084747Z-001 (2)/jNNDFT_publish_no_gui/lda.so";
        CInterface demo = (CInterface) Native.load(libName, CInterface.class);
//1238.0868108401366, 1232.8632901584763, ,4.0071993444865173E8,3.97143579929548E8,-80.93042216577088, -81.18562165757888 ,-161.51939676077606

        //150.10182154640623, 150.59358810567068, ,225683.77966338085,228910.05572378746,-0.0637371869182614, -9.57206045375845

        double gamma_s_[] = {3.97143579929548E8,3.97143579929548E8};

         double rho__[] = {1238.0868108401366, 1232.8632901584763};
        double Fx_[] = {0, 0};
        double dFxdn_[] = {8, 8};
        demo.gga_x(rho__, ao.powdot(gamma_s_, 1.), Fx_, dFxdn_, dero, Fx_.length, 131);
        //  Fx_ = ao.mdot(Fx_, Math.pow(2.0, 1.0 / 3.0));
        //Fx_ = ao.adddot(ao.mdot(Fx, -1), Fx_);
        // Fx_ = ao.mdot(Fx_, 1.0 / Math.pow(2.0, 1.0 / 3.0));
//ao.disp(dFxdn_);
        ao.disp(Fx_);
         ao.disp(rho__);
        ao.disp(ao.divdot(Fx_, rho__));
        //ao.disp(ao.mdot(dero,ao.powdot(gamma_s, 0.5)));

    }

}
