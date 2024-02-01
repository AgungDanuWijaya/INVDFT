package jqc;

import java.util.HashMap;
import function.main_function;
import grid.grid;
import tools.matrix_operation;

/**
 *
 * @author Agung Danu Wijaya
 */
public class exc_functional {

    public double Vxc[][];
    public double Exc;

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
                                v_ggr[i_][j_] += Math.pow(10, -8.0) * r_[3]*dFxdn[p] * datagrid[p_][i_] * (ek[p] / Math.pow(besarr, 2)) * datagrid[p_][j_] * P[i_][j_];
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

    public void GGA(double P[][], grid a, main_function kernel, double datagrid[][], double datagridx[][], double datagridy[][], double datagridz[][],
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

        double rho13[] = kernel.mp.powdot(rho, 1.0 / 3.0);
        double x[] = kernel.mp.divdot(kernel.mp.powdot(gamma_s, 0.5), rho13);
        x = kernel.mp.divdot(x, rho);
        double b = 0.0042;
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
        dFxdn = kernel.mp.mdot(dFxdn, kernel.mp.adddot(b88_g, kernel.mp.mdot(kernel.mp.mdot(x, b88_dg), -1)));

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
        double Fx[] = kernel.ex_ann.Exc(kernel, rho,gamma_s);
        //==============================================		
        double[] dFxdn = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), kernel.ex_ann.Exc(kernel, kernel.mp.adddot(rho, delta),gamma_s)),
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
             double datagridx_[][], double datagridy_[][], double datagridz_[][],HashMap<Integer, HashMap<Integer, double[]>> points) {
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

         double gamax_[] = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.adddot(rhox, kernel.mp.mdot(rho, -2)),rhox_), 1.0 / (kernel.delta_gama*kernel.delta_gama));
        double gamay_[] = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.adddot(rhoy, kernel.mp.mdot(rho, -2)),rhoy_), 1.0 / (kernel.delta_gama*kernel.delta_gama));
        double gamaz_[] = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.adddot(rhoz, kernel.mp.mdot(rho, -2)),rhoz_), 1.0 / (kernel.delta_gama*kernel.delta_gama));
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
        double Fx[] = kernel.ex_ann.Exc_meta(kernel, rho,gamma_s,gamma_s_);
        //==============================================		
        double[] dFxdn = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), kernel.ex_ann.Exc_meta(kernel, kernel.mp.adddot(rho, delta),gamma_s,gamma_s_)),
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
        public void ann_metaGGA_(double P[][], grid a, main_function kernel, double datagrid[][], double datagridx[][], double datagridy[][], double datagridz[][],
         double datagridx_[][], double datagridy_[][], double datagridz_[][],
            HashMap<Integer, HashMap<Integer, double[]>> points) {
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

        double gamax_[] = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.adddot(rhox, kernel.mp.mdot(rho, -2)),rhox_), 1.0 / (kernel.delta_gama*kernel.delta_gama));
        double gamay_[] = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.adddot(rhoy, kernel.mp.mdot(rho, -2)),rhoy_), 1.0 / (kernel.delta_gama*kernel.delta_gama));
        double gamaz_[] = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.adddot(rhoz, kernel.mp.mdot(rho, -2)),rhoz_), 1.0 / (kernel.delta_gama*kernel.delta_gama));
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
        //define your functional========================
        double Fx[] = kernel.ex_ann.Exc_meta(kernel, rho,gamma_s,gamma_s_);
        //==============================================		
        double[] dFxdn = kernel.mp.mdot(kernel.mp.adddot(kernel.mp.mdot(Fx, -1), kernel.ex_ann.Exc(kernel, kernel.mp.adddot(rho, delta),gamma_s)),
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
