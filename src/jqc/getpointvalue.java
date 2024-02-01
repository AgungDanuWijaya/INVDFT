package jqc;

import function.main_function;
import jqc.getdata.datakHF;

/**
 *
 * @author Agung Danu Wijaya
 */
public class getpointvalue {

	main_function kernel;

	public getpointvalue(main_function kernel) {
		this.kernel = kernel;
	}

	public double GTO(double indeks, double r[], double l[]) {
		double besarr = Math.pow(r[0], 2) + Math.pow(r[1], 2) + Math.pow(r[2], 2);
		return Math.pow(r[0], l[0]) * Math.pow(r[1], l[1]) * Math.pow(r[2], l[2]) * Math.exp(-indeks * besarr);
	}

	public double sumbasis(double ca[], double la[], double alphaa[], double Ra[]) {
		double sum = 0;
		for (int m = 0; m < ca.length; m++) {
			double a = alphaa[m];
			double res = GTO(a, Ra, la);
			sum += ca[m] * res * kernel.cn.norm(a, la);
		}
		return sum;
	}

        
        public double op_kin(double ca[], double la[], double alphaa[], double Ra[]) {
        double sum = 0;
        double dx = 0.001;
        double R[] = {Ra[0], Ra[1], Ra[2]};
        for (int m = 0; m < ca.length; m++) {
            double a = alphaa[m];
            double dfx_1 = GTO(a, R, la);
            R[0] += dx;
            double dfx_2 = -2 * GTO(a, R, la);
            R[0] += dx;
            double dfx_3 = GTO(a, R, la);
            double d2fx = dfx_1 + dfx_2 + dfx_3;
//System.out.println( dfx_1+"==="+dfx_2+" "+dfx_3+"  "+d2fx);
            R[0] = Ra[0];
            double dfy_1 = GTO(a, R, la);
            R[1] += dx;
            double dfy_2 = -2 * GTO(a, R, la);
            R[1] += dx;
            double dfy_3 = GTO(a, R, la);
            double d2fy = dfy_1 + dfy_2 + dfy_3;

            R[1] = Ra[1];
            double dfz_1 = GTO(a, R, la);
            R[2] += dx;
            double dfz_2 = -2 * GTO(a, R, la);
            R[2] += dx;
            double dfz_3 = GTO(a, R, la);
            double d2fz = (dfz_1 + dfz_2 + dfz_3);

            double res = d2fx + d2fy + d2fz;
            res=res/(dx*dx);
            
            sum += ca[m] * res * kernel.cn.norm(a, la);
        }
        return sum;
    }

    public double kinetic(datakHF a, double R[]) {
        double r[] = kernel.mp.adddot(R, kernel.mp.mdot(a.R, -1));
        return op_kin(a.c, a.l, a.alpa, r);
    }
	public double bf(datakHF a, double R[]) {
		double r[] = kernel.mp.adddot(R, kernel.mp.mdot(a.R, -1));
		return sumbasis(a.c, a.l, a.alpa, r);
	}
}
