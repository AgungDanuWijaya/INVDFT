package jqc;

import com.sun.jna.Native;
import java.util.HashMap;
import function.main_function;
import grid.grid;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import tools.matrix_operation;

/**
 *
 * @author Agung Danu Wijaya
 */
public class exc_functional {

    public double Vxc[][];
    public double Exc;
    public double gama[][];

    String libName = "/home/agung/project/Quantum-20211024T084747Z-001 (2)/jNNDFT_publish_no_gui/lda.so";

    public void gen(grid a, main_function kernel, double datagrid[][],
            HashMap<Integer, HashMap<Integer, double[]>> points, double Fx[], double dFxdn[], double dFxdgam[], double dFxdngamab[], double datagridx[][], double datagridy[][], double datagridz[][], double gama[][], double gama_b[][]) {

        //double[] dFxdn = kernel.mp.mdot(rho3, (4. / 3.) * fac);
        HashMap<Integer, double[]> pointmap = a.pointsmap(points);
        double Vxc[][] = new double[datagrid[0].length][datagrid[0].length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    double RG[] = pointmap.get(p);
                    Vxc[i][j] += RG[3] * dFxdn[p] * datagrid[p][i] * datagrid[p][j];
                    double amx = (datagridx[p][i] * datagridx[p][j] - datagrid[p][i] * datagrid[p][j]) / kernel.delta_gama;
                    double amy = (datagridy[p][i] * datagridy[p][j] - datagrid[p][i] * datagrid[p][j]) / kernel.delta_gama;
                    double amz = (datagridz[p][i] * datagridz[p][j] - datagrid[p][i] * datagrid[p][j]) / kernel.delta_gama;
                    Vxc[i][j] += 2 * RG[3] * dFxdgam[p] * (gama[0][p] * amx + gama[1][p] * amy + gama[2][p] * amz);
                    Vxc[i][j] += RG[3] * dFxdngamab[p] * (gama_b[0][p] * amx + gama_b[1][p] * amy + gama_b[2][p] * amz);

                }
            }
        }
        double w[] = new double[datagrid.length];

        for (int p = 0; p < datagrid.length; p++) {
            double RG[] = pointmap.get(p);
            w[p] = RG[3];
        }
        double Exc = kernel.mp.sum(kernel.mp.mdot(w, kernel.mp.mdot(Fx, 2)));
        this.Exc = Exc;
        this.Vxc = Vxc;
    }

    public double[][] gam(double P[][], grid a, main_function kernel, double datagrid[][], double datagridx[][], double datagridy[][], double datagridz[][],
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        double rho[] = new double[datagrid.length];
        double rhox[] = new double[datagrid.length];
        double rhoy[] = new double[datagrid.length];
        double rhoz[] = new double[datagrid.length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    rho[p] += datagrid[p][i] * datagrid[p][j] * P[i][j];
                    rhox[p] += datagridx[p][i] * datagridx[p][j] * P[i][j];
                    rhoy[p] += datagridy[p][i] * datagridy[p][j] * P[i][j];
                    rhoz[p] += datagridz[p][i] * datagridz[p][j] * P[i][j];
                }
            }
        }

        double gamax[] = kernel.mp.mdot(kernel.mp.adddot(rhox, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamay[] = kernel.mp.mdot(kernel.mp.adddot(rhoy, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamaz[] = kernel.mp.mdot(kernel.mp.adddot(rhoz, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamma[][] = {gamax, gamay, gamaz};
        double gamma_s[] = kernel.mp.powdot(gamax, 2);
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(gamay, 2));
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(gamaz, 2));
        double[][] gam = new double[3][];
        gam[0] = rho;
        gam[1] = gamma_s;
        this.gama = gamma;

        return gam;
    }

    public void GGA_libxc(int in, double gam[][], double gamc[][], double P[][], grid a, main_function kernel, double datagrid[][], double datagridx[][], double datagridy[][], double datagridz[][],
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        double rho[] = new double[datagrid.length];
        double rhox[] = new double[datagrid.length];
        double rhoy[] = new double[datagrid.length];
        double rhoz[] = new double[datagrid.length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    rho[p] += datagrid[p][i] * datagrid[p][j] * P[i][j];
                    rhox[p] += datagridx[p][i] * datagridx[p][j] * P[i][j];
                    rhoy[p] += datagridy[p][i] * datagridy[p][j] * P[i][j];
                    rhoz[p] += datagridz[p][i] * datagridz[p][j] * P[i][j];
                }
            }
        }

        double gamax[] = kernel.mp.mdot(kernel.mp.adddot(rhox, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamay[] = kernel.mp.mdot(kernel.mp.adddot(rhoy, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamaz[] = kernel.mp.mdot(kernel.mp.adddot(rhoz, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamma[][] = {gamax, gamay, gamaz};
        double gamma_s[] = kernel.mp.powdot(gamax, 2);
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(gamay, 2));
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(gamaz, 2));

        double Fx[] = new double[rho.length];
        // Loading dynamically the library
        CInterface demo = (CInterface) Native.load(libName, CInterface.class);

        double dFxdn[] = new double[Fx.length];
        double dero[] = new double[Fx.length];

        for (int i = 0; i < rho.length; i++) {
            double Fx_[] = {0, 0};
            double dFxdn_[] = {0, 0};
            double rho_[] = {rho[i], 0};
            double dero_[] = {rho[i], 0};
            double gamas_s_[] = {gamma_s[i], 0};
            demo.gga_x(rho_, gamas_s_, Fx_, dFxdn_, dero_, Fx_.length, 106);
            double rho__[] = {gam[0][i], gamc[0][i]};
            Fx[i] = Fx_[0];
            dFxdn[i] = dFxdn_[0];

            double gamas_s__[] = {gam[1][i], gamc[1][i]};
            demo.gga_x(rho__, gamas_s__, Fx_, dFxdn_, dero_, Fx_.length, 131);

            Fx[i] += 0.5 * Fx_[0];

            dFxdn[i] += dFxdn_[in];
        }

        //gen(a, kernel, datagrid, points, Fx, dFxdn);
    }

    public void GGA_lyp(int in, double gama[][], double gama_b[][], double gam[][], double gamc[][], double P[][], grid a, main_function kernel, double datagrid[][], double datagridx[][], double datagridy[][], double datagridz[][],
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        double gamma_s[] = kernel.mp.powdot(kernel.mp.adddot(gama[0], gama_b[0]), 2);
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(kernel.mp.adddot(gama[1], gama_b[1]), 2));
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(kernel.mp.adddot(gama[2], gama_b[2]), 2));

        double Fx[] = new double[gam[0].length];
        // Loading dynamically the library
        CInterface demo = (CInterface) Native.load(libName, CInterface.class);

        double dFxdn[] = new double[Fx.length];
        double dero[] = new double[Fx.length];

        for (int i = 0; i < gam[0].length; i++) {
            double Fx_[] = {0, 0};
            double dFxdn_[] = {0, 0};
            double dero_[] = {gam[0][i], 0};

            double rho__[] = {gam[0][i] + 1E-20, gamc[0][i] + 1E-20};
            double gamas_s__[] = {gam[1][i], gamma_s[i], gamc[1][i]};
            demo.gga_x(rho__, gamas_s__, Fx_, dFxdn_, dero_, Fx_.length, 131);
            double f = (Fx_[0] * (rho__[0] + rho__[1]) / rho__[0]);
            Fx[i] = f;

            dFxdn[i] = dFxdn_[in];
        }

        //gen(a, kernel, datagrid, points, Fx, dFxdn);
    }

    public void GGA_C(int in, double gama[][], double gama_b[][], double gam[][], double gamc[][], double P[][], grid a, main_function kernel, double datagrid[][], double datagridx[][], double datagridy[][], double datagridz[][],
            HashMap<Integer, HashMap<Integer, double[]>> points) {

        double gamma_s[] = kernel.mp.powdot(kernel.mp.adddot(gama[0], gama_b[0]), 2);
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(kernel.mp.adddot(gama[1], gama_b[1]), 2));
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(kernel.mp.adddot(gama[2], gama_b[2]), 2));
        gamc[0] = kernel.mp.adddot(gamc[0], 1E-250);
        gam[0] = kernel.mp.adddot(gam[0], 1E-250);
        gamc[1] = kernel.mp.adddot(gamc[1], 1E-250);
        gam[1] = kernel.mp.adddot(gam[1], 1E-250);
        gamma_s = kernel.mp.adddot(gamma_s, 1E-250);
        //define your functional========================
        double Fx[] = kernel.ex_ann.Exc_cor(in, kernel, gam[0], gamc[0], gam[1], gamc[1], gamma_s);
        //==============================================		
        if (in == 0) {
            double del[] = kernel.mp.mdot(gam[1], 1e-25);
            double[] rt1 = kernel.ex_ann.Exc_cor(in, kernel, gam[0], gamc[0], kernel.mp.adddot(gam[1], del), gamc[1], gamma_s);

            double[] dFxdgam = kernel.mp.divdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), rt1),
                    del);
            double del1[] = kernel.mp.mdot(gamma_s, 1e-25);
            double[] rt2 = kernel.ex_ann.Exc_cor(in, kernel, gam[0], gamc[0], gam[1], gamc[1], kernel.mp.adddot(gamma_s, del1));

            double[] dFxdgam2 = kernel.mp.divdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), rt2),
                    del1);
            double delta[] = kernel.mp.mdot(gam[0], 1e-9);
            double[] dFxdn = kernel.mp.divdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), kernel.ex_ann.Exc_cor(in, kernel, kernel.mp.adddot(gam[0], delta), gamc[0], gam[1], gamc[1], gamma_s)),
                    delta);

            gen(a, kernel, datagrid, points, Fx, dFxdn, dFxdgam, dFxdgam2, datagridx, datagridy, datagridz, gama, gama_b);
        } else {
            double del[] = kernel.mp.mdot(gamc[1], 1e-25);
            double[] rt1 = kernel.ex_ann.Exc_cor(in, kernel, gam[0], gamc[0], gam[1], kernel.mp.adddot(gamc[1], del), gamma_s);

            double[] dFxdgam = kernel.mp.divdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), rt1),
                    del);
            double del1[] = kernel.mp.mdot(gamma_s, 1e-25);
            double[] rt2 = kernel.ex_ann.Exc_cor(in, kernel, gam[0], gamc[0], gam[1], gamc[1], kernel.mp.adddot(gamma_s, del1));

            double[] dFxdgam2 = kernel.mp.divdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), rt2),
                    del1);

            double delta[] = kernel.mp.mdot(gamc[0], 1e-9);
            double[] dFxdn = kernel.mp.divdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), kernel.ex_ann.Exc_cor(in, kernel, gam[0], kernel.mp.adddot(gamc[0], delta), gam[1], gamc[1], gamma_s)),
                    delta);

            gen(a, kernel, datagrid, points, Fx, dFxdn, dFxdgam, dFxdgam2, datagridx, datagridy, datagridz, gama, gama_b);
        }
    }

    public void LDA(double P[][], grid a, main_function kernel, double datagrid[][],
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        double rho[] = new double[datagrid.length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    rho[p] += datagrid[p][i] * datagrid[p][j] * P[i][j];
                }
            }
        }

        double alpha = 2.0 / 3.0;
        double fac = -2.25 * alpha * Math.pow(0.75 / Math.PI, 1. / 3.);
        double rho3[] = kernel.mp.powdot(rho, 1. / 3.);
        double Fx[] = kernel.mp.mdot(kernel.mp.mdot(rho, rho3), fac);
        double[] dFxdn = kernel.mp.mdot(rho3, (4. / 3.) * fac);

        HashMap<Integer, double[]> pointmap = a.pointsmap(points);
        double Vxc[][] = new double[datagrid[0].length][datagrid[0].length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    double RG[] = pointmap.get(p);
                    Vxc[i][j] += RG[3] * dFxdn[p] * datagrid[p][i] * datagrid[p][j];
                }
            }
        }
        double w[] = new double[datagrid.length];
        for (int p = 0; p < datagrid.length; p++) {
            double RG[] = pointmap.get(p);
            w[p] = RG[3];
        }
        double Exc = kernel.mp.sum(kernel.mp.mdot(w, kernel.mp.mdot(Fx, 2)));
        this.Exc = Exc;
        this.Vxc = Vxc;
    }

    public void LDA_PZ(double P[][], grid a, main_function kernel, double datagrid[][],
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        double rho[] = new double[datagrid.length];
        double exc_result[][] = new double[2][datagrid.length];
        double cor_result[][] = new double[2][datagrid.length];
        exchange exc = new exchange();
        correlation cor = new correlation();
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    rho[p] += datagrid[p][i] * datagrid[p][j] * P[i][j];

                }
            }
            if (rho[p] > Math.pow(10, -11)) {
                double rs = Math.pow(3.0 / (4.0 * Math.PI), 1.0 / 3.0) / Math.pow(rho[p], 1.0 / 3.0);
                double ex[] = exc.slater(rs, rho[p]);
                double ec[] = cor.pz(rs, rho[p], 1);
                //  exc_result[0][p] = ex[0];
                // exc_result[1][p] = ex[1];
                cor_result[0][p] = ec[0];
                cor_result[1][p] = ec[1];
                //   System.out.println("  exc_result[0][p]"+  cor_result[0][p]+" "+  cor_result[1][p]);
            }
        }

        double alpha = 2.0 / 3.0;
        double fac = -2.25 * alpha * Math.pow(0.75 / Math.PI, 1. / 3.);
        double rho3[] = kernel.mp.powdot(rho, 1. / 3.);
        double Fx_[] = kernel.mp.mdot(kernel.mp.mdot(rho, rho3), fac);
        double[] dFxdn_ = kernel.mp.mdot(rho3, (4. / 3.) * fac);

        double Fx[] = kernel.mp.adddot(Fx_, cor_result[0]);
        double[] dFxdn = kernel.mp.adddot(dFxdn_, cor_result[1]);

        HashMap<Integer, double[]> pointmap = a.pointsmap(points);
        double Vxc[][] = new double[datagrid[0].length][datagrid[0].length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    double RG[] = pointmap.get(p);
                    Vxc[i][j] += RG[3] * dFxdn[p] * datagrid[p][i] * datagrid[p][j];
                }
            }
        }
        double w[] = new double[datagrid.length];
        for (int p = 0; p < datagrid.length; p++) {
            double RG[] = pointmap.get(p);
            w[p] = RG[3];
        }
        double Exc = kernel.mp.sum(kernel.mp.mdot(w, kernel.mp.mdot(Fx, 2)));
        this.Exc = Exc;
        this.Vxc = Vxc;

    }

    public void GGR1(double P[][], grid a, main_function kernel, double datagrid[][], double datagrid_kinetik[][],
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        double rho[] = new double[datagrid.length];

        HashMap<Integer, double[]> pointmap = a.pointsmap(points);
        double ek[] = new double[datagrid.length];

        for (int p = 0; p < datagrid.length; p++) {
            double r[] = pointmap.get(p);
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    ek[p] += r[3] * datagrid[p][i] * datagrid_kinetik[p][j] * P[i][j];
                    rho[p] += datagrid[p][i] * datagrid[p][j] * P[i][j];
                    //   System.out.println(datagrid_kinetik[p][j]);
                }
            }
        }
        double alpha = 2.0 / 3.0;
        double fac = -2.25 * alpha * Math.pow(0.75 / Math.PI, 1. / 3.);
        double rho3[] = kernel.mp.powdot(rho, 1. / 3.);
        double Fx[] = kernel.mp.mdot(kernel.mp.mdot(rho, rho3), fac);
        double[] dFxdn = kernel.mp.mdot(rho3, (4. / 3.) * fac);

        double etot = 0;
        double v_ggr[][] = new double[datagrid[0].length][datagrid[0].length];
        for (int p_ = 0; p_ < datagrid.length; p_++) {
            double r_[] = pointmap.get(p_);
            for (int p = 0; p < datagrid.length; p++) {
                double r[] = pointmap.get(p);
                double besarr = Math.pow(r_[0] - r[0], 2) + Math.pow(r_[1] - r[1], 2) + Math.pow(r_[2] - r[2], 2);
                besarr = Math.pow(besarr, 0.5);

                // if (besarr > Math.pow(1.0/Math.abs(rho[p]), 1.0/3.0)) {
                if (besarr > 0) {
                    if (besarr < 1) {
                        for (int i_ = 0; i_ < datagrid[p_].length; i_++) {
                            for (int j_ = 0; j_ < datagrid[p_].length; j_++) {
                                v_ggr[i_][j_] += Math.pow(10, -8.0) * r_[3] * dFxdn[p] * datagrid[p_][i_] * (ek[p] / Math.pow(besarr, 2)) * datagrid[p_][j_] * P[i_][j_];
                                etot += v_ggr[i_][j_];
                            }
                        }
                    }
                }

            }
        }
        // this.Exc = etot;
        //this.Vxc = v_ggr;
        this.Exc += etot;
        this.Vxc = new matrix_operation().adddot(this.Vxc, v_ggr);

    }

    public void ann(double P[][], grid a, main_function kernel, double datagrid[][],
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        double rho[] = new double[datagrid.length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    rho[p] += datagrid[p][i] * datagrid[p][j] * P[i][j];
                }
            }
        }

        double delta = Math.pow(10, -7);
        //define your functional========================
        double Fx[] = kernel.ex_ann.Ex(kernel, rho);
        //==============================================		
        double[] dFxdn = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), kernel.ex_ann.Ex(kernel, kernel.mp.adddot(rho, delta))),
                1.0 / delta);

        HashMap<Integer, double[]> pointmap = a.pointsmap(points);
        double Vxc[][] = new double[datagrid[0].length][datagrid[0].length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    double RG[] = pointmap.get(p);
                    Vxc[i][j] += RG[3] * dFxdn[p] * datagrid[p][i] * datagrid[p][j];
                }
            }
        }
        double w[] = new double[datagrid.length];
        for (int p = 0; p < datagrid.length; p++) {
            double RG[] = pointmap.get(p);
            w[p] = RG[3];
        }
        double Exc = kernel.mp.sum(kernel.mp.mdot(w, kernel.mp.mdot(Fx, 2)));
        this.Exc = Exc;
        this.Vxc = Vxc;
    }

    double asinh(double x) {
        return Math.log(x + Math.sqrt(x * x + 1.0));
    }

    public double[][] b88(double rho[], double gamma_s[], main_function kernel) {
        double rho13[] = kernel.mp.powdot(rho, 1.0 / 3.0);

        double x[] = kernel.mp.divdot(kernel.mp.powdot(gamma_s, 0.5), rho13);
        x = kernel.mp.divdot(x, rho);
        double b = 0.0042;
        //double b = kernel.weight[0][0][0][0];
        double b88_g[] = new double[x.length];
        double b88_dg[] = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            double b88_g_ = -1.5 * Math.pow(3.0 / 4.0 / Math.PI, 1.0 / 3.0) - b * x[i] * x[i] / (1.0 + 6.0 * b * x[i] * asinh(x[i]));
            double num = (6.0 * b * b * x[i] * x[i] * (x[i] / Math.sqrt(x[i] * x[i] + 1.0) - asinh(x[i]))) - (2.0 * b * x[i]);
            double denom = Math.pow(1.0 + (6.0 * b * x[i] * asinh(x[i])), 2.0);
            double b88_dg_ = num / denom;
            b88_g[i] = b88_g_;
            b88_dg[i] = b88_dg_;
        }

        double Fx[] = kernel.mp.mdot(kernel.mp.mdot(rho13, rho), b88_g);

        double[] dFxdn = kernel.mp.mdot(rho13, (4.0 / 3.0));
        double[] dFxdgam = kernel.mp.divdot(kernel.mp.mdot(b88_dg, (0.5)), kernel.mp.powdot(gamma_s, 0.5));
        dFxdn = kernel.mp.mdot(dFxdn, kernel.mp.adddot(b88_g, kernel.mp.mdot(kernel.mp.mdot(x, b88_dg), -1)));
        double tr[][] = {Fx, dFxdn, dFxdgam};
        return tr;
    }

    public double[][] b88_(double rho[], double gamma_s[], main_function kernel) {
        double rho13[] = kernel.mp.powdot(rho, 1.0 / 3.0);

        double x[] = kernel.mp.divdot(kernel.mp.powdot(gamma_s, 0.5), rho13);
        x = kernel.mp.divdot(x, rho);
        double b = 0.0042;
        //double b = kernel.weight[0][0][0][0];
        double b88_g[] = new double[x.length];
        double b88_dg[] = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            double b88_g_ = -1.5 * Math.pow(3.0 / 4.0 / Math.PI, 1.0 / 3.0) - b * x[i] * x[i] / (1.0 + 6.0 * b * x[i] * asinh(x[i]));
            double num = (6.0 * b * b * x[i] * x[i] * (x[i] / Math.sqrt(x[i] * x[i] + 1.0) - asinh(x[i]))) - (2.0 * b * x[i]);
            double denom = Math.pow(1.0 + (6.0 * b * x[i] * asinh(x[i])), 2.0);
            double b88_dg_ = num / denom;
            b88_g[i] = b88_g_;
            b88_dg[i] = b88_dg_;
        }

        double Fx[] = kernel.mp.mdot(kernel.mp.mdot(rho13, rho), b88_g);

        double[] dFxdn = kernel.mp.mdot(rho13, (4.0 / 3.0));
        double[] dFxdgam = kernel.mp.divdot(kernel.mp.mdot(b88_dg, (0.5)), kernel.mp.powdot(gamma_s, 0.5));
        dFxdn = kernel.mp.mdot(dFxdn, kernel.mp.adddot(b88_g, kernel.mp.mdot(kernel.mp.mdot(x, b88_dg), -1)));
        double tr[][] = {Fx, dFxdn, dFxdgam};
        return tr;
    }


    
    public void GGA_X(String xc, double P[][], grid a, main_function kernel, double datagrid[][], double datagridx[][], double datagridy[][], double datagridz[][],
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        double rho[] = new double[datagrid.length];
        double rhox[] = new double[datagrid.length];
        double rhoy[] = new double[datagrid.length];
        double rhoz[] = new double[datagrid.length];

        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    rho[p] += datagrid[p][i] * datagrid[p][j] * P[i][j];
                    rhox[p] += datagridx[p][i] * datagridx[p][j] * P[i][j];
                    rhoy[p] += datagridy[p][i] * datagridy[p][j] * P[i][j];
                    rhoz[p] += datagridz[p][i] * datagridz[p][j] * P[i][j];

                }
            }
        }

        double gamax[] = kernel.mp.mdot(kernel.mp.adddot(rhox, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamay[] = kernel.mp.mdot(kernel.mp.adddot(rhoy, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamaz[] = kernel.mp.mdot(kernel.mp.adddot(rhoz, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamma[][] = {gamax, gamay, gamaz};

        double gamma_s[] = kernel.mp.powdot(gamax, 2);
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(gamay, 2));
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(gamaz, 2));

        double Fx[] = null;
        double[] dFxdn = null;
        double[] dFxdgam = null;

        if (xc.equals("B88")) {
            double[][] rt = b88(rho, gamma_s, kernel);
            Fx = rt[0];
            dFxdn = rt[1];
            double del[] = kernel.mp.mdot(gamma_s, 1e-5);
            double[][] rt1 = b88(rho, kernel.mp.adddot(gamma_s, del), kernel);

            dFxdgam = kernel.mp.divdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), rt1[0]),
                    del);
            
            
       
        
    
        } else if (xc.equals("ann_gga")) {
            rho = kernel.mp.adddot(rho, 1E-250);
            gamma_s = kernel.mp.adddot(gamma_s, 1E-250);
            double delta[] = kernel.mp.mdot(rho, 1e-9);

            //define your functional========================
            Fx = kernel.ex_ann.Exc(kernel, rho, gamma_s);
            //==============================================		
            dFxdn = kernel.mp.divdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), kernel.ex_ann.Exc(kernel, kernel.mp.adddot(rho, delta), gamma_s)),
                    delta);
            double del[] = kernel.mp.mdot(gamma_s, 1e-25);
            double Fx1[] = kernel.ex_ann.Exc(kernel, rho, kernel.mp.adddot(gamma_s, del));
            dFxdgam = kernel.mp.divdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), Fx1),
                    del);
            
          ObjectOutputStream outputStream;
      /*  try {
            outputStream = new ObjectOutputStream(new FileOutputStream(kernel.URL_ANN+"_rho"));
            outputStream.writeObject(rho);
        } catch (FileNotFoundException ex) {
            System.err.println(ex);
        } catch (IOException ex) {
            System.err.println(ex);
        }
         try {
            outputStream = new ObjectOutputStream(new FileOutputStream(kernel.URL_ANN+"_gamma_s"));
            outputStream.writeObject(gamma_s);
        } catch (FileNotFoundException ex) {
            System.err.println(ex);
        } catch (IOException ex) {
            System.err.println(ex);
        }*/
         
        }

        HashMap<Integer, double[]> pointmap = a.pointsmap(points);
        double Vxc[][] = new double[datagrid[0].length][datagrid[0].length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    double RG[] = pointmap.get(p);

                    Vxc[i][j] += RG[3] * dFxdn[p] * datagrid[p][i] * datagrid[p][j];

                    double amx = (datagridx[p][i] * datagridx[p][j] - datagrid[p][i] * datagrid[p][j]) / kernel.delta_gama;
                    double amy = (datagridy[p][i] * datagridy[p][j] - datagrid[p][i] * datagrid[p][j]) / kernel.delta_gama;
                    double amz = (datagridz[p][i] * datagridz[p][j] - datagrid[p][i] * datagrid[p][j]) / kernel.delta_gama;
                    Vxc[i][j] += 2 * RG[3] * dFxdgam[p] * (gamax[p] * amx + gamay[p] * amy + gamaz[p] * amz);

                    // System.out.println(dFxdgam[p]+" "+dFxdn[p]+" "+(2 * RG[3] * dFxdgam[p] * taugamma_s[p]* datagrid[p][i] * datagrid[p][j])+" "+RG[3] * dFxdn[p] * datagrid[p][i] * datagrid[p][j]);
                }
            }

        }
        // System.out.println(sum);

        double w[] = new double[datagrid.length];
        for (int p = 0; p < datagrid.length; p++) {
            double RG[] = pointmap.get(p);
            w[p] = RG[3];
        }
        double Exc = kernel.mp.sum(kernel.mp.mdot(w, kernel.mp.mdot(Fx, 2)));
        this.Exc = Exc;
        this.Vxc = Vxc;

    }

    public void ann_GGA_1(double P[][], grid a, main_function kernel, double datagrid[][], double datagridx[][], double datagridy[][], double datagridz[][],
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        double rho[] = new double[datagrid.length];
        double rhox[] = new double[datagrid.length];
        double rhoy[] = new double[datagrid.length];
        double rhoz[] = new double[datagrid.length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    rho[p] += datagrid[p][i] * datagrid[p][j] * P[i][j];
                    rhox[p] += datagridx[p][i] * datagridx[p][j] * P[i][j];
                    rhoy[p] += datagridy[p][i] * datagridy[p][j] * P[i][j];
                    rhoz[p] += datagridz[p][i] * datagridz[p][j] * P[i][j];
                }
            }
        }

        double gamax[] = kernel.mp.mdot(kernel.mp.adddot(rhox, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamay[] = kernel.mp.mdot(kernel.mp.adddot(rhoy, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamaz[] = kernel.mp.mdot(kernel.mp.adddot(rhoz, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamma[][] = {gamax, gamay, gamaz};
        double gamma_s[] = kernel.mp.powdot(gamax, 2);
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(gamay, 2));
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(gamaz, 2));

        double delta = Math.pow(10, -7);
        //define your functional========================
        double Fx[] = kernel.ex_ann.Exc__(kernel, rho, gamma_s);
        //==============================================		
        double[] dFxdn = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), kernel.ex_ann.Exc__(kernel, kernel.mp.adddot(rho, delta), gamma_s)),
                1.0 / delta);

        HashMap<Integer, double[]> pointmap = a.pointsmap(points);
        double Vxc[][] = new double[datagrid[0].length][datagrid[0].length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    double RG[] = pointmap.get(p);
                    Vxc[i][j] += RG[3] * dFxdn[p] * datagrid[p][i] * datagrid[p][j];
                }
            }
        }
        double w[] = new double[datagrid.length];
        for (int p = 0; p < datagrid.length; p++) {
            double RG[] = pointmap.get(p);
            w[p] = RG[3];
        }
        double Exc = kernel.mp.sum(kernel.mp.mdot(w, kernel.mp.mdot(Fx, 2)));
        this.Exc = Exc;
        this.Vxc = Vxc;

    }

    public void ann_GGA(double P[][], grid a, main_function kernel, double datagrid[][], double datagridx[][], double datagridy[][], double datagridz[][],
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        double rho[] = new double[datagrid.length];
        double rhox[] = new double[datagrid.length];
        double rhoy[] = new double[datagrid.length];
        double rhoz[] = new double[datagrid.length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    rho[p] += datagrid[p][i] * datagrid[p][j] * P[i][j];
                    rhox[p] += datagridx[p][i] * datagridx[p][j] * P[i][j];
                    rhoy[p] += datagridy[p][i] * datagridy[p][j] * P[i][j];
                    rhoz[p] += datagridz[p][i] * datagridz[p][j] * P[i][j];
                }
            }
        }

        double gamax[] = kernel.mp.mdot(kernel.mp.adddot(rhox, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamay[] = kernel.mp.mdot(kernel.mp.adddot(rhoy, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamaz[] = kernel.mp.mdot(kernel.mp.adddot(rhoz, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamma[][] = {gamax, gamay, gamaz};
        double gamma_s[] = kernel.mp.powdot(gamax, 2);
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(gamay, 2));
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(gamaz, 2));

        double delta = Math.pow(10, -7);
        //define your functional========================
        double Fx[] = kernel.ex_ann.Exc(kernel, rho, gamma_s);
        //==============================================		
        double[] dFxdn = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), kernel.ex_ann.Exc(kernel, kernel.mp.adddot(rho, delta), gamma_s)),
                1.0 / delta);

        HashMap<Integer, double[]> pointmap = a.pointsmap(points);
        double Vxc[][] = new double[datagrid[0].length][datagrid[0].length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    double RG[] = pointmap.get(p);
                    Vxc[i][j] += RG[3] * dFxdn[p] * datagrid[p][i] * datagrid[p][j];
                }
            }
        }
        double w[] = new double[datagrid.length];
        for (int p = 0; p < datagrid.length; p++) {
            double RG[] = pointmap.get(p);
            w[p] = RG[3];
        }
        double Exc = kernel.mp.sum(kernel.mp.mdot(w, kernel.mp.mdot(Fx, 2)));
        this.Exc = Exc;
        this.Vxc = Vxc;

    }

    public void ann_metaGGA(double P[][], grid a, main_function kernel, double datagrid[][], double datagridx[][], double datagridy[][], double datagridz[][],
            double datagridx_[][], double datagridy_[][], double datagridz_[][], HashMap<Integer, HashMap<Integer, double[]>> points) {
        double rho[] = new double[datagrid.length];
        double rhox[] = new double[datagrid.length];
        double rhoy[] = new double[datagrid.length];
        double rhoz[] = new double[datagrid.length];

        double rhox_[] = new double[datagrid.length];
        double rhoy_[] = new double[datagrid.length];
        double rhoz_[] = new double[datagrid.length];

        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    rho[p] += datagrid[p][i] * datagrid[p][j] * P[i][j];
                    rhox[p] += datagridx[p][i] * datagridx[p][j] * P[i][j];
                    rhoy[p] += datagridy[p][i] * datagridy[p][j] * P[i][j];
                    rhoz[p] += datagridz[p][i] * datagridz[p][j] * P[i][j];

                    rhox_[p] += datagridx_[p][i] * datagridx_[p][j] * P[i][j];
                    rhoy_[p] += datagridy_[p][i] * datagridy_[p][j] * P[i][j];
                    rhoz_[p] += datagridz_[p][i] * datagridz_[p][j] * P[i][j];
                }
            }
        }

        double gamax_[] = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.adddot(rhox, kernel.mp.mdot(rho, -2)), rhox_), 1.0 / (kernel.delta_gama * kernel.delta_gama));
        double gamay_[] = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.adddot(rhoy, kernel.mp.mdot(rho, -2)), rhoy_), 1.0 / (kernel.delta_gama * kernel.delta_gama));
        double gamaz_[] = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.adddot(rhoz, kernel.mp.mdot(rho, -2)), rhoz_), 1.0 / (kernel.delta_gama * kernel.delta_gama));
        double gamma_[][] = {gamax_, gamay_, gamaz_};
        double gamma_s_[] = kernel.mp.powdot(gamax_, 2);
        gamma_s_ = kernel.mp.adddot(gamma_s_, kernel.mp.powdot(gamay_, 2));
        gamma_s_ = kernel.mp.adddot(gamma_s_, kernel.mp.powdot(gamaz_, 2));

        double gamax[] = kernel.mp.mdot(kernel.mp.adddot(rhox, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamay[] = kernel.mp.mdot(kernel.mp.adddot(rhoy, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamaz[] = kernel.mp.mdot(kernel.mp.adddot(rhoz, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamma[][] = {gamax, gamay, gamaz};
        double gamma_s[] = kernel.mp.powdot(gamax, 2);
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(gamay, 2));
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(gamaz, 2));

        double delta = Math.pow(10, -7);
        //  System.out.println("jqc.exc_functional.ann_metaGGA()");
        //define your functional========================
        double Fx[] = kernel.ex_ann.Exc_meta(kernel, rho, gamma_s, gamma_s_);
        //==============================================		
        double[] dFxdn = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), kernel.ex_ann.Exc_meta(kernel, kernel.mp.adddot(rho, delta), gamma_s, gamma_s_)),
                1.0 / delta);

        HashMap<Integer, double[]> pointmap = a.pointsmap(points);
        double Vxc[][] = new double[datagrid[0].length][datagrid[0].length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    double RG[] = pointmap.get(p);
                    Vxc[i][j] += RG[3] * dFxdn[p] * datagrid[p][i] * datagrid[p][j];
                }
            }
        }
        double w[] = new double[datagrid.length];
        for (int p = 0; p < datagrid.length; p++) {
            double RG[] = pointmap.get(p);
            w[p] = RG[3];
        }
        double Exc = kernel.mp.sum(kernel.mp.mdot(w, kernel.mp.mdot(Fx, 2)));
        this.Exc = Exc;
        this.Vxc = Vxc;

    }

    public void ann_GGA_(double P[][], grid a, main_function kernel, double datagrid[][], double datagridx[][], double datagridy[][], double datagridz[][],
            double datagridx_[][], double datagridy_[][], double datagridz_[][], HashMap<Integer, HashMap<Integer, double[]>> points) {
        double rho[] = new double[datagrid.length];
        double rhox[] = new double[datagrid.length];
        double rhoy[] = new double[datagrid.length];
        double rhoz[] = new double[datagrid.length];

        double rhox_[] = new double[datagrid.length];
        double rhoy_[] = new double[datagrid.length];
        double rhoz_[] = new double[datagrid.length];

        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    rho[p] += datagrid[p][i] * datagrid[p][j] * P[i][j];
                    rhox[p] += datagridx[p][i] * datagridx[p][j] * P[i][j];
                    rhoy[p] += datagridy[p][i] * datagridy[p][j] * P[i][j];
                    rhoz[p] += datagridz[p][i] * datagridz[p][j] * P[i][j];

                    rhox_[p] += datagridx_[p][i] * datagridx_[p][j] * P[i][j];
                    rhoy_[p] += datagridy_[p][i] * datagridy_[p][j] * P[i][j];
                    rhoz_[p] += datagridz_[p][i] * datagridz_[p][j] * P[i][j];
                }
            }
        }

        double gamax_[] = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.adddot(rhox, kernel.mp.mdot(rho, -2)), rhox_), 1.0 / (kernel.delta_gama * kernel.delta_gama));
        double gamay_[] = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.adddot(rhoy, kernel.mp.mdot(rho, -2)), rhoy_), 1.0 / (kernel.delta_gama * kernel.delta_gama));
        double gamaz_[] = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.adddot(rhoz, kernel.mp.mdot(rho, -2)), rhoz_), 1.0 / (kernel.delta_gama * kernel.delta_gama));
        double gamma_[][] = {gamax_, gamay_, gamaz_};
        double gamma_s_[] = kernel.mp.powdot(gamax_, 2);
        gamma_s_ = kernel.mp.adddot(gamma_s_, kernel.mp.powdot(gamay_, 2));
        gamma_s_ = kernel.mp.adddot(gamma_s_, kernel.mp.powdot(gamaz_, 2));

        double gamax[] = kernel.mp.mdot(kernel.mp.adddot(rhox, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamay[] = kernel.mp.mdot(kernel.mp.adddot(rhoy, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamaz[] = kernel.mp.mdot(kernel.mp.adddot(rhoz, kernel.mp.mdot(rho, -1)), 1.0 / kernel.delta_gama);
        double gamma[][] = {gamax, gamay, gamaz};
        double gamma_s[] = kernel.mp.powdot(gamax, 2);
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(gamay, 2));
        gamma_s = kernel.mp.adddot(gamma_s, kernel.mp.powdot(gamaz, 2));

        double delta = Math.pow(10, -7);
        //  System.out.println("jqc.exc_functional.ann_metaGGA()");
        //define your functional========================
        double Fx[] = kernel.ex_ann.Exc_meta(kernel, rho, gamma_s, kernel.mp.mdot(gamma_s_, 0));
        //==============================================		
        double[] dFxdn = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), kernel.ex_ann.Exc_meta(kernel, kernel.mp.adddot(rho, delta), gamma_s, gamma_s_)),
                1.0 / delta);

        HashMap<Integer, double[]> pointmap = a.pointsmap(points);
        double Vxc[][] = new double[datagrid[0].length][datagrid[0].length];
        for (int p = 0; p < datagrid.length; p++) {
            for (int i = 0; i < datagrid[p].length; i++) {
                for (int j = 0; j < datagrid[p].length; j++) {
                    double RG[] = pointmap.get(p);
                    Vxc[i][j] += RG[3] * dFxdn[p] * datagrid[p][i] * datagrid[p][j];
                }
            }
        }
        double w[] = new double[datagrid.length];
        for (int p = 0; p < datagrid.length; p++) {
            double RG[] = pointmap.get(p);
            w[p] = RG[3];
        }
        double Exc = kernel.mp.sum(kernel.mp.mdot(w, kernel.mp.mdot(Fx, 2)));
        this.Exc = Exc;
        this.Vxc = Vxc;

    }

}