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
public class Interface {

    public String method = "nr";
    public String thr_nri = "1000";
    public String thr_eri = "80000";
    public String out = "";
    public String url_db = "jdbc:mysql://127.0.0.1:3306/";
    public String user_db = "mabok_janda";
    public String pass_db = "yut28092018DAM^";
    public double proses_int = 0;
    public int print = 0;
    public String ann_conf = "1,1";
    public String data_src = "MySQL";//[MySQL,Manual]
    public String dens = "2";
    public String basis = "6-31g";
    public boolean int_stat = false;
    public String conv = "0.000001";
    public String mix = "0.75";
    public String exc_tipe = "ANN";
    public String base = "/home/agung/project/Quantum-20211024T084747Z-001 (2)/Quantum (copy)/JQC_data";
    public String name_exc = "ann_new_dft";
    public String tugas = "1";
    public String st = "0.1";
    public int cluster_num = 1;
    public String data_geo = "{\n" +
"\"name\": [\"H\",\"OH,O,H}\",\"H2,H,H}\",\"Na2,Na,Na}\",\"H2S,S,H,H}\",\"H2O2,H,H,O,O}\"],\n" +
"\"con\": 0.0015936254980079682,\n" +
"\"re\": [\n" +
"	[-0.5],[-0.162231075697211],[-0.164621513944223], [-0.0264541832669323],[-0.27601593625498],[-0.40207171314741]\n" +
"                          ],\n" +
"\"simpul\": [1,1],\n" +
"\"name_ann\": \"ann_new_dft\"\n" +
"}";
    public double pa_ne[] = {1.0, 0.000000001};
    public String user_param = "1";

    public static void main(String[] args) {
        Interface in = new Interface();
        try {
            ann_dft_hf a = new ann_dft_hf(in);
            String method = "DFT";
            if (method.equals("DFT_HF")) {
                a.energi();
            } else if (method.equals("DFT")) {
                double r = a.energi_dft();
            } else if (method.equals("UHF")) {
                a.energi_hf();
            }
        } catch (Exception e) {
            System.out.println(e);
        }
    }
}
