/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package jqc;

import com.sun.jna.Native;
import tools.matrix_operation;

/**
 *
 * @author agung
 */
public class test_lyp {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        double A = 0.04918;
        double B = 0.132;
        double C = 0.2533;
        double Dd = 0.349;

        double CF = (3.0 / 10.0) * Math.pow(3.0 * Math.pow(Math.PI, 2), 2.0 / 3.0);

        double gamma_s_[] = {1.3221842650466303E-2,1.4334852236073482E-2};
//80.937395877003 0.0
//2.6001395205033247E-13 2.706573597620902E-13 1.3221842650466303E-24 1.4334852236073482E-24 5.509090235749842E-24 -Infinity
        double rho__[] = {0.1653514803498147,0.1653514803498147};
        double a = rho__[0];
        double b = rho__[1];
        double gaa_ = gamma_s_[0];
        double gbb_ = gamma_s_[1];

        double n = a + b;
        double gnn_ =5.509090235749842E-2; 
        double gaa = Math.pow(gaa_, 1);
        double gbb = Math.pow(gbb_, 1);
        double gnn = Math.pow(gnn_, 1);

        double icbrtn = Math.pow(n, -1.0 / 3.0);
        double P = 1.0 / (1.0 + Dd * icbrtn);
        double omega = Math.exp(-C * icbrtn) * P * Math.pow(n, -11.0 / 3.0);
        double delta = icbrtn * (C + Dd * P);
        System.out.println(omega+" "+icbrtn);
        double n2 = n * n;
        double f = -A * (4.0 * a * b * P / n
                + B * omega
                * (a * b * (Math.pow(2.0, 11.0 / 3.0) * CF
                * (Math.pow(a, 8.0 / 3.0) + Math.pow(b, 8.0 / 3.0))
                + (47.0 - 7.0 * delta) * gnn / 18.0
                - (2.5 - delta / 18.0) * (gaa + gbb)
                - ((delta - 11.0) / 9.0) * (a * gaa + b * gbb) / n)
                - (2.0 / 3.0) * n2 * gnn + (2.0 / 3.0 * n2 - a * a) * gbb
                + (2.0 / 3.0 * n2 - b * b) * gaa));

        System.out.println(f * rho__[0] / (a + b) + " " + (f / (a + b))+" "+f);
        
        String libName = "/home/agung/project/Quantum-20211024T084747Z-001 (2)/jNNDFT_publish_no_gui/lda.so";
        CInterface demo = (CInterface) Native.load(libName, CInterface.class);
        matrix_operation ao = new matrix_operation();

        double Fx_[] = {0, 0};

        double dero[] = {0, 0};
        double dFxdn_[] = {0, 0};
          //    double gamma_s__[] = {75,5*4*3,4*4*3};
double gamma_s[] = {1.3221842650466303E-2,5.509090235749842E-2,1.4334852236073482E-2};
        
        demo.gga_x(rho__, ao.powdot(gamma_s, 1.), Fx_, dFxdn_, dero, Fx_.length, 131);
        
        
        ao.disp(Fx_);
        
        System.out.println(Fx_[0]*(a + b)/rho__[0]);
      
        ao.disp(ao.divdot(Fx_, rho__));
     //   System.err.println(Double.isInfinite(Double.in));
    }

}
