package ann;

import ann.ann;
import function.main_function;

public class drv_ann {

    public double[] Ex(main_function kernel, double input[]) {
        double Ex_ann[] = new double[input.length];
        ann c = new ann();
        for (int i = 0; i < input.length; i++) {
            double[] rho = {input[i]};
            Ex_ann[i] = kernel.weight[2][0][0][1] * Math.pow(input[i], kernel.weight[2][0][0][0]);
            //Ex_ann[i] = c.f_ANN(rho, kernel.c_node, kernel.weight, kernel.thv, kernel.th);
        }
        return Ex_ann;
    }

    public double[] Exc(main_function kernel, double input[], double gama[]) {
        double Ex_ann[] = new double[input.length];
        ann c = new ann();
        for (int i = 0; i < input.length; i++) {
            double[] rho = {Math.pow(input[i], kernel.weight[2][0][0][0]) * Math.pow(gama[i], kernel.weight[2][0][0][1])};

            if (input[i] < Math.pow(10, -11)) {
                double[] o = {0};

                rho = o;
            } else if (gama[i] < Math.pow(10, -11)) {
                double[] o = {0};

                rho = o;
            } else if (input[i] < Math.pow(10, -11) & gama[i] < Math.pow(10, -11)) {
                double[] o = {0};

                rho = o;
            }

            Ex_ann[i] = c.f_ANN(rho, kernel.c_node, kernel.weight, kernel.thv, kernel.th);

        }
        return Ex_ann;
    }

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

    public double[] Exc_meta_(main_function kernel, double input[], double gama[], double gamas[]) {
        double Ex_ann[] = new double[input.length];
        ann c = new ann();
        for (int i = 0; i < input.length; i++) {
            double[] rho = {Math.pow(input[i], kernel.weight[2][0][0][0]) * Math.pow(gama[i], kernel.weight[2][0][0][1]), Math.pow(input[i], kernel.weight[2][0][0][2]) * Math.pow(gamas[i], kernel.weight[2][0][0][3])};

            if (input[i] < Math.pow(10, -11)) {
                double[] o = {0, 0};

                rho = o;
            } else if (gama[i] < Math.pow(10, -11)) {
                double[] o = {0, 0};

                rho = o;
            } else if (input[i] < Math.pow(10, -11) & gama[i] < Math.pow(10, -11)) {
                double[] o = {0, 0};

                rho = o;
            }

            Ex_ann[i] = c.f_ANN(rho, kernel.c_node, kernel.weight, kernel.thv, kernel.th);

        }
        return Ex_ann;
    }
}
