# INVDFT: Java Program to optimize parameter for Exchange Correlation Functional in  Density Functional Theory
INVDFT, a Java based software, employs the SCF technique to solve the KSDFT equation. The optimization of XC function parameters in INVDFT is accomplished using the Newton-Raphson and Monte Carlo methods
## Features
- Parameter optimization of Exchange Correlation Functional
- SCF with electron density mixing
- DFT: Neural Network, LDA and GGA functionals
- Hartree-Fock
## Requirements
- Linux OS
- Java (https://www.java.com/en/download/) version >= 11
-  MySQL
## Example of parameter optimization

$\epsilon_{xc}=\beta \int(\rho_{\-}^\gamma+\rho_+^\gamma) dr^3$

------------
	public double[] Ex(main_function kernel, double input[]) {
	        double Ex_ann[] = new double[input.length];
	        ann c = new ann();
	        for (int i = 0; i < input.length; i++) {
	            double[] rho = {input[i]};
	            Ex_ann[i] = kernel.weight[2][0][0][1] * Math.pow(input[i], kernel.weight[2][0][0][0]);
	        }
	        return Ex_ann;
	}
------------
------------
	{
	"name": ["H","OH,O,H}","H2,H,H}","Na2,Na,Na}","H2S,S,H,H}"],
	"con": 0.0015936254980079682,
	"re": [
	[-0.5],[-0.162231075697211],[-0.164621513944223],[-0.0264541832669323],[-0.27601593625498]
                          ],
	"simpul": [1,1],
	"name_ann": "ann_new_dft"
	}
------------

------------
    public String method = "nr";
    public String thr_nri = "1000";
    public String thr_eri = "80000";
    public String out = "";
    public String url_db = "jdbc:mysql://127.0.0.1:3306/";
    public String user_db = "user_db";
    public String pass_db = "pass_db";
    public double proses_int = 0;
    public int print = 0;
    public String ann_conf = "1,1";
    public String data_src = "MySQL";//[MySQL,Manual]
    public String dens = "0";
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
    public String data_geo = "{\n"
            + "\"name\": [\"H\",\"OH,O,H}\",\"H2,H,H}\",\"Na2,Na,Na}\",\"H2S,S,H,H}\",\"H2O2,H,H,O,O}\"],\n"
            + "\"con\": 0.0015936254980079682,\n"
            + "\"re\": [\n"
            + "	[-0.5],[-0.162231075697211],[-0.164621513944223], [-0.0264541832669323],[-0.27601593625498],[-0.40207171314741]\n"
            + "                          ],\n"
            + "\"simpul\": [1,1],\n"
            + "\"name_ann\": \"ann_new_dft\"\n"
            + "}";
    public double pa_ne[] = {1.0, 0.000000001};
    public String user_param = "1";
    
------------

## Example of DFT Calculation
------------
	H2O;
	Water;
	{"xyz": 
	{"1": [0.0000, 0.0000, 0.1173], 
	"2": [0.0000, 0.7572, -0.4692],
	"3": [0.0000, -0.7572, -0.4692]}, 
	"atom": ["O", "H","H"],
	"Spin_dn": 5, "Spin_up": 5};
------------
