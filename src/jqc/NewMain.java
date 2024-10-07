/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package jqc;

import function.main_function;

/**
 *
 * @author agung
 */
public class NewMain {

    public static double f(double x, double y) {
        return x * x +  x * y  -3;
    }

    public static void main(String[] args) {
        double x = 0;
        double y = 0;
        for (int i = 0; i < 1000; i++) {
            double derx = (f(x + 0.00001, y) - f(x, y)) / 0.00001;
            double dery = (f(x, y + 0.00001) - f(x, y)) / 0.00001;
            x = x - 0.2 * derx;
            y = y - 0.2 * dery;
            System.out.println(i + " : "+x+" "+y+ " " + f(x, y));
        }
    }

}
