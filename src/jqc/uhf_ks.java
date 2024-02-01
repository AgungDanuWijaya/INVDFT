package jqc;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.Map;
import Jama.Matrix;
import function.main_function;
import geo_molecule.data_geo;
import grid.grid;
import java.util.HashMap;

/**
 *
 * @author Agung Danu Wijaya, Dedy Farhamsa
 */
public class uhf_ks {

    main_function kernel;

    public uhf_ks(main_function kernel) {
        this.kernel = kernel;
    }

    public double SCF(String geo) throws IOException, ClassNotFoundException, InterruptedException {
        Map<String, data_geo> data;
        data = kernel.geo.data;
        grid a = new grid(geo, kernel);
        HashMap<Integer, HashMap<Integer, double[]>> points = a.points;
        Map<Integer, getdata.datakHF> bfs = kernel.gdata.get(geo);
        double[][] datagrid = a.setbfamps(kernel, bfs, points);
        double[][] datagrid_kinetik = a.set_gr(kernel, bfs, points);
        double[][] datagrid_x = a.setbfamps_x(kernel, bfs, points);
        double[][] datagrid_y = a.setbfamps_y(kernel, bfs, points);
        double[][] datagrid_z = a.setbfamps_z(kernel, bfs, points);
        double[][] datagrid_x_ = a.setbfamps_x_(kernel, bfs, points);
        double[][] datagrid_y_ = a.setbfamps_y_(kernel, bfs, points);
        double[][] datagrid_z_ = a.setbfamps_z_(kernel, bfs, points);
          double Ej = 0, Exc = 0, Eh = 0;
        double Ejc = 0, Excc = 0, Ehc = 0;
        double S[][] = null;
        double T[][] = null;
        double V[][] = null;
        double G[][][][] = null;
        if (kernel.status_int == 1) {
            kernel.intg.one(geo);
            kernel.intg.two(geo);
            G = kernel.intg.ints;
            S = kernel.intg.S;
            T = kernel.intg.EK;
            V = kernel.intg.EV;
        } else {
            try {
                ObjectInputStream inputStream = new ObjectInputStream(new FileInputStream(kernel.URL_int + "_G"));
                G = (double[][][][]) inputStream.readObject();
                inputStream.close();
                inputStream = new ObjectInputStream(new FileInputStream(kernel.URL_int + "_S"));
                S = (double[][]) inputStream.readObject();
                inputStream.close();
                inputStream = new ObjectInputStream(new FileInputStream(kernel.URL_int + "_T"));
                T = (double[][]) inputStream.readObject();
                inputStream.close();
                inputStream = new ObjectInputStream(new FileInputStream(kernel.URL_int + "_V"));
                V = (double[][]) inputStream.readObject();
                inputStream.close();

            } catch (FileNotFoundException ex) {
                System.out.println(ex);
            } catch (IOException ex) {
                System.out.println(ex);
            }
        }
        double H[][] = kernel.mp.adddot(V, T);
        double C[][] = kernel.gev.gev_run(S, H);
        double kali = kernel.c_mix;
        double En = 0;
        int panjangu = 0;
        int panjangd = 0;
        panjangu = kernel.spinup;
        panjangd = kernel.spindn;
        double DU[][] = new double[C.length][panjangu];
        double DD[][] = new double[C.length][panjangd];
        for (int i = 0; i < C.length; i++) {
            for (int j = 0; j < panjangu; j++) {
                DU[i][j] = C[i][j] + Math.random() * 0.00000001;
            }
            for (int j = 0; j < panjangd; j++) {
                DD[i][j] = C[i][j] + Math.random() * 0.00000001;
            }
        }

        Matrix UU = new Matrix(DU);
        Matrix DUB = UU.times(UU.transpose());
        Matrix UD = new Matrix(DD);
        Matrix DDB = UD.times(UD.transpose());
        double PU[][] = DUB.getArray();
        double PD[][] = DDB.getArray();
        double CUold[][] = kernel.mp.copy(PU);
        double CDold[][] = kernel.mp.copy(PD);
        double enold = 0;
        int stop = 0;
        while (stop == 0) {
            double RPU[][][] = new double[2][S.length][S.length];
            double RPD[][][] = new double[2][S.length][S.length];
            for (int i = 0; i < S.length; i++) {
                for (int j = 0; j < S.length; j++) {
                    for (int k = 0; k < S.length; k++) {
                        for (int l = 0; l < S.length; l++) {
                            RPU[0][i][j] += G[i][j][k][l] * PU[k][l];
                            RPU[1][i][j] += G[i][k][j][l] * PU[k][l];
                            RPD[0][i][j] += G[i][j][k][l] * PD[k][l];
                            RPD[1][i][j] += G[i][k][j][l] * PD[k][l];
                        }
                    }
                }
            }
            
             exc_functional FX = new exc_functional();
            exc_functional FXc = new exc_functional();
            if (kernel.Ex.equals("LDA")) {
                FX.LDA(PU, a, kernel, datagrid, points);
                FXc.LDA(PD, a, kernel, datagrid, points);
            } else if (kernel.Ex.equals("ANN")) {
                FX.ann(PU, a, kernel, datagrid, points);
                FXc.ann(PD, a, kernel, datagrid, points);
            } else if (kernel.Ex.equals("LDA_PZ")) {
                FX.LDA_PZ(PU, a, kernel, datagrid, points);
                FXc.LDA_PZ(PD, a, kernel, datagrid, points);
            }  else if (kernel.Ex.equals("GGA_X_B88")) {
                FX.GGA(PU, a, kernel, datagrid, datagrid_x, datagrid_y, datagrid_z, points);
                FXc.GGA(PD, a, kernel, datagrid, datagrid_x, datagrid_y, datagrid_z, points);
            }else if (kernel.Ex.equals("ann_GGA")) {
                FX.ann_GGA(PU, a, kernel, datagrid, datagrid_x, datagrid_y, datagrid_z, points);
                FXc.ann_GGA(PD, a, kernel, datagrid, datagrid_x, datagrid_y, datagrid_z, points);
            }else if (kernel.Ex.equals("ann_metaGGA")) {
                FX.ann_metaGGA(PU, a, kernel, datagrid, datagrid_x, datagrid_y, datagrid_z,datagrid_x_, datagrid_y_, datagrid_z_, points);
                FXc.ann_metaGGA(PD, a, kernel, datagrid, datagrid_x, datagrid_y, datagrid_z,datagrid_x_, datagrid_y_, datagrid_z_, points);
            }
            

            double Vxc[][] = FX.Vxc;
            Exc = FX.Exc;
            double Vxcc[][] = FXc.Vxc;
            Excc = FXc.Exc;

            double mix=0.0;
            double FU[][] = kernel.mp.adddot(H,
                    kernel.mp.adddot(RPU[0], kernel.mp.adddot(kernel.mp.mdot(RPU[1], -mix), RPD[0])));
            double FD[][] = kernel.mp.adddot(H,
                    kernel.mp.adddot(RPD[0], kernel.mp.adddot(kernel.mp.mdot(RPD[1], -mix), RPU[0])));
            
             double FU_[][] = kernel.mp.adddot(H,
                    kernel.mp.adddot(RPU[0], kernel.mp.adddot(kernel.mp.mdot(RPU[1], -mix), RPD[0])));
            double FD_[][] = kernel.mp.adddot(H,
                    kernel.mp.adddot(RPD[0], kernel.mp.adddot(kernel.mp.mdot(RPD[1], -mix), RPU[0])));
        FU =  kernel.mp.adddot(FU, Vxc);
        FD =  kernel.mp.adddot(FD, Vxcc);
            
            
            
            double CU[][] = kernel.gev.gev_run(S, FU);
            double CD[][] = kernel.gev.gev_run(S, FD);
            for (int i = 0; i < CU.length; i++) {
                for (int j = 0; j < panjangu; j++) {
                    DU[i][j] = CU[i][j];
                }
            }
            for (int i = 0; i < CD.length; i++) {
                for (int j = 0; j < panjangd; j++) {
                    DD[i][j] = CD[i][j];
                }
            }

            UU = new Matrix(DU);
            DUB = UU.times(UU.transpose());
            UD = new Matrix(DD);
            DDB = UD.times(UD.transpose());
            PU = DUB.getArray();
            PD = DDB.getArray();

            PU = kernel.mp.adddot(kernel.mp.mdot(PU, 1 - kali), kernel.mp.mdot(CUold, kali));
            PD = kernel.mp.adddot(kernel.mp.mdot(PD, 1 - kali), kernel.mp.mdot(CDold, kali));
            CUold = kernel.mp.copy(PU);
            CDold = kernel.mp.copy(PD);
            En = kernel.mp.sum(kernel.mp.mdot(kernel.mp.adddot(H, FU_), PU));
            En += kernel.mp.sum(kernel.mp.mdot(kernel.mp.adddot(H, FD_), PD));
            En +=Exc;
            En+=Excc;
            En /= 2.0;
            if (Math.abs(En - enold) < Double.parseDouble(kernel.a.conv)) {
                stop = 1;
            } else {
                enold = En;
            }
            System.out.println(String.valueOf(En + kernel.geo.energi(data.get(geo))) + "\n");

        }
        return En + kernel.geo.energi(data.get(geo));
    }

}
