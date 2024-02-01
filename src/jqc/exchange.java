/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jqc;

/**
 *
 * @author agung
 */
public class exchange {

    public double[] slater(double rs, double rho) {
        //double alpha = 2.0 / 3.0;
        // double fac = -2.25 * alpha * Math.pow(0.75 / Math.PI, 1. / 3.);
        double f = -9.0 / 8.0 * Math.pow(3.0 / (2.0 * Math.PI), 2.0 / 3.0);
        double alpha = 2.0 / 3.0;
        double ex = rho * f * alpha / rs;
        double vx = (4.0 / 3.0) * f * (alpha / rs);
        double result[] = {ex, vx};
        return result;
    }

}
