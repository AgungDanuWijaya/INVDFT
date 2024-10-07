/* Source : Computational Techniques in Quantum Chemistry
 * and Molecular Physics
 * see equation number 35 on page 366
*/
package integral;

import function.main_function;

/**
 *
 * @author Agung Danu Wijaya
 */
public class f {

	main_function a;

	public f(main_function kernel) {
		this.a = kernel;
	}

	public double lgamma(double x) {
		x = a.specialFunc.gamma(x);
		return Math.log(Math.abs(x));
	}

	public double[] _gser(double a, double x) {
		double ITMAX = 100;
		double EPS = 3.0E-7;

		double gln = lgamma(a);
		if (x == 0) {
			double res[] = { 0, gln };
			return res;
		}
		double ap = a;
		double delt = 1.0 / a;
		double sum = delt;
		for (int i = 0; i < ITMAX; i++) {
			ap = ap + 1;
			delt = delt * x / ap;
			sum = sum + delt;
			if (Math.abs(delt) < Math.abs(sum) * EPS) {
				i = (int) (ITMAX + 90);
			}
		}
		double gamser = sum * Math.exp(-x + a * Math.log(x) - gln);
		double res[] = { gamser, gln };
		return res;
	}

	public double[] _gcf(double a, double x) {
		double ITMAX = 100;
		double EPS = 3.0E-7;
		double FPMIN = 1.0E-30;

		double gln = lgamma(a);
		double b = x + 1.0 - a;
		double c = 1.0 / FPMIN;
		double d = 1.0 / b;
		double h = d;
		for (int i = 1; i < ITMAX + 1; i++) {
			double an = -i * (i - a);
			b += 2;
			d = an * d + b;
			if (Math.abs(d) < FPMIN) {
				d = FPMIN;
			}
			c = b + an / c;
			if (Math.abs(c) < FPMIN) {
				c = FPMIN;
			}
			d = 1 / d;
			double delt = d * c;
			h *= delt;
			if (Math.abs(delt - 1) < EPS) {
				i = (int) (ITMAX + 90);
			}
		}
		double gammcf = Math.exp(-x + a * Math.log(x) - gln) * h;

		double res[] = { gammcf, gln };
		return res;
	}

	public double gamm_inc(double a, double x) {
		double res[];
		if (x < a + 1) {
			res = _gser(a, x);
		} else {
			res = _gcf(a, x);
			res[0] = 1 - res[0];
		}
		return Math.exp(res[1]) * res[0];
	}

	public double clc_F(double m, double x) {
		double SMALL = 1E-12;
		x = Math.max(x, SMALL);
		return 0.5 * Math.pow(x, -m - 0.5) * gamm_inc(m + 0.5, x);
	}

	
}