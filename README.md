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
## Instalation 
-  Mysql configuration
  - Instal Mysql, for detail step you can refer tp this website https://www.digitalocean.com/community/tutorials/how-to-install-mysql-on-ubuntu-20-04
  - Import Dump20240206.sql file to database
-  Netbeans configuration
  - Install netbeans from https://netbeans.apache.org/front/main/download/index.html
  - Import project INVDFT to netbeans project
## Example of LSDA XC parameter optimization
Misalkan kita mempunyai model LSDA XC sebagai berikut : $\epsilon_{xc}=\beta \int(\rho_{\-}^\gamma+\rho_+^\gamma) dr^3$. Kita  akan mencari $\beta$ dan $\gamma$ dengan menggunkaan INVDFT. Langkah pertama adalah mendefinisikan model XC agar dikenali oleh INVDFT. Definisi model xc dapat dilakukan pada drv_ann.java pada ann package.
Untuk model model LSDA, potongan script akhir dari drv_ann.java adalah sebagai berikut.

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
Setelah kita mendefinisikan model XC, kita akan memberikan data ajar yang digunakan untuk fitting parameter dalam model LSDA. Penginputan data ajar dapat dilakukan pada variabel data_geo dalam file Interface.java di Interface package. Data ajr yang digunakan terdiri dari Atomisasi energi dan energi dari molekul. Pada contoh dibawah kita menggunakan "H"=-0.5 Hartree sebagau data ajar enrgi molekul dan "OH,O,H}"=-0.162231075697211 Hartree sebagai data ajar atomisasi energi dari molekul OH.

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

Setalah melakukan penginputan data ajar pada variabel data_geo, script final dari file Interface.java menjadi berikut.

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
Setelah semua persiapan selesai dilakukan, kita dapat memulai proses optimisasi. Langkah-langkah dari proses optimisasi adalah sebagai berikut.
- Run file inversi.java on Interface package
- Run file cluster.java on Interface package
- run script SELECT * FROM Quantum.error on mysql server to see error on every iteration; 

## Example of Neural-Like XC parameter optimization
Misal kita akan membuat model XC dari gambar dibawah.

<img src="https://github.com/AgungDanuWijaya/INVDFT/blob/main/Screenshot%20from%202024-02-06%2015-20-15.png" alt="dftk logo" height="200px" />

Langkah pertama adalah mendefinisikan model neural like XC agar dikenali oleh INVDFT. Definisi model xc dapat dilakukan pada drv_ann.java pada ann package. Untuk model model neural like, potongan script akhir dari drv_ann.java adalah sebagai berikut. Pada script dibawah input berperan sebagai $\rho$, gamma berperan sebagai $\nabla \rho$, dan gammas berperan sebagai $\nabla^2 \rho$.



------------
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
------------

Kemudian menginput data ajar pada file interface.java. Data ajar diinput pada variabel data_geo. Contoh script dapat dilihat dibawah.

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
    public String ann_conf = "3,2,2,1";
    public String data_src = "MySQL";//[MySQL,Manual]
    public String dens = "0";
    public String basis = "6-31g";
    public boolean int_stat = false;
    public String conv = "0.000001";
    public String mix = "0.75";
    public String exc_tipe = "ann_metaGGA";//ann_metaGGA//ANN
    public String base = "/home/agung/project/Quantum-20211024T084747Z-001 (2)/jNNDFT_publish_no_gui/JQC_data";
    public String name_exc = "ann_new_dft";
    public String tugas = "1";
    public String st = "0.1";
    public int cluster_num = 5;
    public String data_geo = "{\n"
            + "\"name\": [\"H\",\"OH,O,H}\",\"H2,H,H}\",\"Na2,Na,Na}\",\"H2S,S,H,H}\",\"H2O2,H,H,O,O}\"],\n"
            + "\"con\": 0.0015936254980079682,\n"
            + "\"re\": [\n"
            + "	[-0.5],[-0.162231075697211],[-0.164621513944223], [-0.0264541832669323],[-0.27601593625498],[-0.40207171314741]\n"
            + "                          ],\n"
            + "\"simpul\": [3,2,2,1],\n"
            + "\"name_ann\": \"ann_new_dft\"\n"
            + "}";
    public double pa_ne[] = {1.0, 0.000000001,1.0, 0.000000001};
    public String user_param = "1";
   
    
------------

Langkah berikutnya sama dengan langkah pada contoh LSDA XC parameter optimization

## Example of DFT Calculation

Edit data_geo on file Interface.java in Interface package

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

run file run_dft.java in Interface package

## Save Integral

Edit data_geo on file Interface.java in Interface package

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

run file save_int.java in Interface package

