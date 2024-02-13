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
  - Install Mysql, for detail step you can refer tp this website https://www.digitalocean.com/community/tutorials/how-to-install-mysql-on-ubuntu-20-04
  - Import Dump20240206.sql file to database
-  Netbeans configuration
  - Install netbeans from https://netbeans.apache.org/front/main/download/index.html
  - Import project INVDFT to netbeans project
## Example of LSDA XC parameter optimization
Suppose we have the LSDA XC model as follows: $\epsilon_{xc}=\beta \int(\rho_{-}^\gamma+\rho_+^\gamma) dr^3$. We will find $\beta$ and $\gamma$ using INVDFT. The first step is to define the XC model to be recognized by INVDFT. The definition of the xc model can be done in drv_ann.java in the ann package.
For the LSDA model, the final script section of drv_ann.java is as follows.

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
The initialization of kernel.weight[2][0][0][0] and kernel.weight[2][0][0][1] in the example above is done in the variable pa_ne[] = {1.0, 0.000000001} in the Interface.java. So in this example kernel.weight[2][0][0][0]=1.0 and kernel.weight[2][0][0][1]=0.000000001.
After defining the XC model, we will provide training data used to fit parameters in the LSDA model. Inputting training data can be done in the variable data_geo in the Interface.java file in the Interface package. The training data used consists of Atomization energy and energy of molecules. In the example below, we use "H"=-0.5 Hartree as the molecular energy training data and "OH,O,H}"=-0.162231075697211 Hartree as the atomization energy training data for the OH molecule. Please remember,  that integrals for all molecular geometries used as teaching data must be saved first. This aims to ensure that integral calculations are not carried out repeatedly in the XC parameter optimization process.

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

After inputting the training data into the data_geo variable, the final script of the Interface.java file becomes as follows:

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

After all preparations have been completed, we can begin the optimization process. The steps of the optimization process are as follows:

- Run the inversi.java file in the Interface package.
- Run the cluster.java file in the Interface package.
- Run the script SELECT * FROM Quantum.error on the MySQL server to see errors on every iteration.

## Example of Neural-Like XC parameter optimization
Suppose we want to create an XC model based on the image below.

<img src="https://github.com/AgungDanuWijaya/INVDFT/blob/main/Screenshot%20from%202024-02-06%2015-20-15.png" alt="dftk logo" height="200px" />

The first step is to define the neural-like XC model to be recognized by INVDFT. The definition of the XC model can be done in drv_ann.java in the ann package. For the neural-like model, the final script section of drv_ann.java is as follows. In the script below, the input acts as $\rho$, gamma acts as $\nabla \rho$, and gammas act as $\nabla^2 \rho$.



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

Then input the training data in the file interface.java. The training data is input into the data_geo variable. An example script can be seen below.

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

The next steps are the same as the steps in the example of LSDA XC parameter optimization.

## Example of DFT Calculation

Edit data_geo , proses_int, and data_src on file Interface.java in Interface package

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

------------
proses_int = 1;
data_src = "Manual";

------------
The variable process_int=1 aims to calculate the integral again according to the geometry in data_geo. Meanwhile, data_src="Manual", aims to ensure that the geometry used in DFT calculations uses the geometry in the data_geo variable, not from geometric data stored in the database.<p>

run file run_dft.java in Interface package<p>
Users can change the XC functional in the exc_tipe variable, where the available functions are LSDA, GGA_B_88, ANN, ann_metaGGA. ANN and ann_metaGGA are user-defined XC types.

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


For more detailed information, you can contact us at email: wijayadanuagung@gmail.com
