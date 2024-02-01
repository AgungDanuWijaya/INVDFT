package jqc;

import java.util.HashMap;
import java.util.Map;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import Jama.Matrix;
import function.main_function;
import geo_molecule.data_geo;
import grid.grid;

/**
 *
 * @author Agung Danu Wijaya, Dedy Farhamsa
 */
public class uks {

    public double SCF(String geo, main_function kernel) throws InterruptedException, ClassNotFoundException {
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


        Map<String, data_geo> data = kernel.geo.data;
        double Ej = 0, Exc = 0, Eh = 0;
        double Ejc = 0, Excc = 0, Ehc = 0;
        double G[][][][] = null;
        double S[][] = null;
        double T[][] = null;
        double V[][] = null;
        int id = 0;
        double totalma[] = new double[2];
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
        double Cc[][] = kernel.gev.gev_run(S, H);
        double Cold[][];
        double Coldc[][];
        double mix_c = kernel.c_mix;
        int panjangc = kernel.spinup;
        int panjango = kernel.spindn;
        double[][] Do = new double[C.length][panjango];
        double[][] Dc = new double[C.length][panjangc];
        randy c=new randy();
        c.init();
        for (int i = 0; i < C.length; i++) {
            for (int j = 0; j < panjango; j++) {
                Do[i][j] = C[i][j] + c.rand()* 0.01;
            }
            for (int j = 0; j < panjangc; j++) {
                Dc[i][j] = Cc[i][j] + c.rand() * 0.01;
            }
        }
        Matrix Uo = new Matrix(Do);
        Matrix DBo = Uo.times(Uo.transpose());
        Matrix Uc = new Matrix(Dc);
        Matrix DBc = Uc.times(Uc.transpose());
        double P[][] = DBo.getArray();
        double Pc[][] = DBc.getArray();
        double enold = 0;
        int stop = 0;
        int itemax = 0;
        while (stop == 0) {
            if (itemax++ > 100) {
             //            System.out.println(Eh +","+ (0.5 * Exc) +","+ Ej+","+Ehc +","+ (0.5 * Excc) +","+ Ejc +","+ kernel.geo.energi(data.get(geo))+",no");
  
                stop = 1;
            }
            double J[][] = new double[S.length][S[0].length];
            double Jc[][] = new double[S.length][S[0].length];
            for (int i = 0; i < S.length; i++) {
                for (int j = 0; j < S.length; j++) {
                    for (int k = 0; k < S.length; k++) {
                        for (int l = 0; l < S.length; l++) {
                            J[i][j] += G[i][j][k][l] * (P[k][l] + Pc[k][l]) * 0.5;
                            Jc[i][j] += G[i][j][k][l] * (P[k][l] + Pc[k][l]) * 0.5;
                        }
                    }
                }
            }

            Ej = kernel.mp.sum(kernel.mp.mdot(J, P));
            Eh = kernel.mp.sum(kernel.mp.mdot(H, P));
            double Ek = kernel.mp.sum(kernel.mp.mdot(T, P));
            double Ev = kernel.mp.sum(kernel.mp.mdot(V, P));
            double total = Eh + 0.5 * Exc + Ej;

            Ejc = kernel.mp.sum(kernel.mp.mdot(Jc, Pc));
            Ehc = kernel.mp.sum(kernel.mp.mdot(H, Pc));

            double totalc = Ehc + 0.5 * Excc + Ejc + kernel.geo.energi(data.get(geo));
            total += totalc;

            totalma[id] = total;
          
           /* if ( kernel . spindn + kernel . spinup == 1) {
       
total = Ek + Ev ;
stop = 1;
totalma [0] = total ;
totalma [1] = total ;
}*/
            id++;
            if (id == 2) {
                id = 0;
            }

            if (Math.abs(total - enold) < Double.parseDouble(kernel.a.conv)) {
                stop = 1;
         //       System.out.println(Eh +","+ (0.5 * Exc) +","+ Ej+","+Ehc +","+ (0.5 * Excc) +","+ Ejc +","+ kernel.geo.energi(data.get(geo)));
            } else {
                enold = total;
            }
            exc_functional FX = new exc_functional();
            exc_functional FXc = new exc_functional();
            if (kernel.Ex.equals("LDA")) {
                FX.LDA(P, a, kernel, datagrid, points);
                FXc.LDA(Pc, a, kernel, datagrid, points);
            } else if (kernel.Ex.equals("ANN")) {
                FX.ann(P, a, kernel, datagrid, points);
                FXc.ann(Pc, a, kernel, datagrid, points);
            } else if (kernel.Ex.equals("LDA_PZ")) {
                FX.LDA_PZ(P, a, kernel, datagrid, points);
                FXc.LDA_PZ(Pc, a, kernel, datagrid, points);
            }  else if (kernel.Ex.equals("GGA_X_B88")) {
                FX.GGA(P, a, kernel, datagrid, datagrid_x, datagrid_y, datagrid_z, points);
                FXc.GGA(Pc, a, kernel, datagrid, datagrid_x, datagrid_y, datagrid_z, points);
            }else if (kernel.Ex.equals("ann_GGA")) {
                FX.ann_GGA(P, a, kernel, datagrid, datagrid_x, datagrid_y, datagrid_z, points);
                FXc.ann_GGA(Pc, a, kernel, datagrid, datagrid_x, datagrid_y, datagrid_z, points);
            }else if (kernel.Ex.equals("ann_metaGGA")) {
                FX.ann_metaGGA(P, a, kernel, datagrid, datagrid_x, datagrid_y, datagrid_z,datagrid_x_, datagrid_y_, datagrid_z_, points);
                FXc.ann_metaGGA(Pc, a, kernel, datagrid, datagrid_x, datagrid_y, datagrid_z,datagrid_x_, datagrid_y_, datagrid_z_, points);
            }
            
            else if (kernel.Ex.equals("GGR")) {
                FX.ann(P, a, kernel, datagrid, points);
                FXc.ann(Pc, a, kernel, datagrid, points);
                if (kernel.spinup + kernel.spindn > 1) {
                    FX.GGR1(P, a, kernel, datagrid, datagrid_kinetik, points);
                    FXc.GGR1(Pc, a, kernel, datagrid, datagrid_kinetik, points);
                }
                //FX.cLDA(P, a, kernel, datagrid, points);
                //FXc.cLDA(Pc, a, kernel, datagrid, points);
            }

            double Vxc[][] = FX.Vxc;
            Exc = FX.Exc;
            double Vxcc[][] = FXc.Vxc;
            Excc = FXc.Exc;

            double RP[][] = kernel.mp.adddot(J, kernel.mp.adddot(Jc, Vxc));
            double F[][] = kernel.mp.adddot(H, RP);
            C = kernel.gev.gev_run(S, F);

            double RPc[][] = kernel.mp.adddot(J, kernel.mp.adddot(Jc, Vxcc));
            double Fc[][] = kernel.mp.adddot(H, RPc);
            Cc = kernel.gev.gev_run(S, Fc);

            for (int i = 0; i < C.length; i++) {
                for (int j = 0; j < panjango; j++) {
                    Do[i][j] = C[i][j];
                }
                for (int j = 0; j < panjangc; j++) {
                    Dc[i][j] = Cc[i][j];
                }
            }
            Uo = new Matrix(Do);
            DBo = Uo.times(Uo.transpose());
            Uc = new Matrix(Dc);
            DBc = Uc.times(Uc.transpose());

            Cold = kernel.mp.copy(P);
            P = DBo.getArray();
            P = kernel.mp.adddot(kernel.mp.mdot(P, 1 - mix_c), kernel.mp.mdot(Cold, mix_c));
            Coldc = kernel.mp.copy(Pc);
            Pc = DBc.getArray();
            Pc = kernel.mp.adddot(kernel.mp.mdot(Pc, 1 - mix_c), kernel.mp.mdot(Coldc, mix_c));
            if (kernel.a.print == 1) {
                System.out.println(String.valueOf(total) + "");
            }

        }
        double Ek = kernel.mp.sum(kernel.mp.mdot(T, P));
        double Ev = kernel.mp.sum(kernel.mp.mdot(V, P));
        double Nai = kernel.geo.energi(data.get(geo));
        double rt[] = {Math.abs(Exc / Math.min(totalma[0], totalma[1])),
            Math.abs(Ej / Math.min(totalma[0], totalma[1])), Math.abs(Ek / Math.min(totalma[0], totalma[1])),
            Math.abs(Ev / Math.min(totalma[0], totalma[1])), Math.abs(Nai / Math.min(totalma[0], totalma[1]))};
        kernel.input_mixer = rt;
        return Math.min(totalma[0], totalma[1]);
    }
}
