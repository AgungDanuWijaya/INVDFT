/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package run;

/**
 *
 * @author agung
 */
public class monte_carlo {

    public static double y(double x[]) {
        double y_ = (x[0] * x[0] - 2 * x[0] + 1)-(x[1] * x[1] - 4 * x[1] + 4);
        return y_;
    }

    public static void main(String args[]) {

        double best = 9999999;
        double x_new[] = {9,9};
        double del = 9;

        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 80; j++) {
                double in[] = new double[x_new.length];
                for (int k = 0; k < x_new.length; k++) {
                    double x = (x_new[k] - del) + Math.random() * (x_new[k] + del - (x_new[k] - del));
                
                    
                    in[k] = x;
                }

                double y_ = y(in);
               // System.out.println(y_);
                if (Math.abs(y_) < best) {
                    best = Math.abs(y_);
                    for (int k = 0; k < in.length; k++) {
                        x_new[k] = in[k];
                    }

                }
            }

            del /= 1.3;

            System.out.println(x_new[0]+" "+x_new[1] + " aa " + y(x_new));
        }

    }
}
