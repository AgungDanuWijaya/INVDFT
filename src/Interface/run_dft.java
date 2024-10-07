/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package Interface;

import run.ann_dft_hf;

/**
 *
 * @author agung
 */
public class run_dft {



    public static void main(String[] args) {
        Interface in = new Interface();
        try {
            ann_dft_hf a = new ann_dft_hf(in);
            String method = "DFT";
            if (method.equals("DFT_HF")) {
               a.energi();
            } else if (method.equals("DFT")) {
                double r = a.energi_dft();
                System.out.println(r);
            } else if (method.equals("UHF")) {
                a.energi_hf();
            }
        } catch (Exception e) {
            System.out.println(e);
        }
    }
}
