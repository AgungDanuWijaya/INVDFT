/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package run;

/**
 *
 * @author agung
 */
public class monte_carlo1 {

    public static double y(double x[]) {
        double y_ = x[0] * x[0] - 4 * x[0] + 4.0;
        y_ /= 2.0 * x[0] * x[0] - 2.0 * x[0] + 1.0;
        return y_ - 1.0 / 13.0;
    }

    public static void main(String args[]) {

        double y = 80;
        double best = 9999999;
        double x_new = 9;
        double del = 1;

        for (int i = 0; i < 80; i++) {
            double x_1 = x_new - del;
            double x_2 = x_new + del;
            for (int j = 0; j < y; j++) {

                double x = x_1 + Math.random() * (x_2 - x_1);

                double in[] = {x};
                double y_ = y(in);
                if (Math.abs(y_) < best) {
                    best = Math.abs(y_);
                    x_new = x;
                }
            }
            // best = 9999999;
            del /= 1.01;
            double x[] = {x_new};
            System.out.println(x_new + " aa " + y(x));
        }

    }
}
