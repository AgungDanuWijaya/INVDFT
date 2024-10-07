package run;

import java.text.SimpleDateFormat;
import java.util.Date;

import function.main_function;

/**
 *
 * @author Agung Danu Wijaya
 */
public class test_integral {

    main_function master = new main_function();

    public void rapat() {
        double a1 = 1.1695961;
        double a2 = 0.380389;
        int l = 0;
        int m = 0;
        int n = 1;
        double l1[] = {l, m, n};
        double l2[] = {l, m, n};
        double ra[] = {0., 0., 0.0916858};
        double rb[] = {0., 0., 0.0916858};
        double r = (master.rhoe.S(a1, a2, l1, l2, ra, rb));
        System.out.println("S " + r);
    }

    public void EK() {
        double a1 = 1.1695961;
        double a2 = 1.1695961;
        int l = 1;
        int m = 1;
        int n = 0;
        double l1[] = {l, m, n};
        double l2[] = {l, m, n};
        double ra[] = {0., 0., 0.04851804};
        double rb[] = {0., 0., 0.04851804};
        System.out.println("EK " + master.T.EK(a1, a2, l1, l2, ra, rb));
    }

    public void NAI() {
        double a1 = 0.38389;
        double a2 = 0.38389;
        int l = 0;
        int m = 0;
        int n = 0;
        double l1[] = {l, m, n};
        double l2[] = {l, m, n};
        double ra[] = {0., 0, 0};
        double rb[] = {0, 0, 0};
        double rc[] = rb;
        double r = master.nai.NAI(a1, a2, l1, l2, ra, rb, rc);
        // System.err.println(r);
        // r=a.naiN.NAIMC(a1, a2, l1, l2, ra, rb, rc);
        double r1 = master.naiN.NAI(a1, a2, l1, l2, ra, rb, rc);
        System.out.println("NAI A " + r + " NAI N " + r1);
    }

    public void ERI() {
        double a1 = 0.130;
        double a2 = 0.130;
        double a3 = 0.23;
        double a4 = 0.23;
        int l = 1;
        int m = 1;
        int n = 1;
        double l1[] = {l, m, n};
        double l2[] = {l, m, n};
        double l3[] = {l, m, n};
        double l4[] = {l, m, n};
        double ra[] = {0.0, 0.0, 2.16858};
        double rb[] = {0.0, 0.0, 2.1858};
        double rc[] = {0.0, 0.0, 0.016858};
        double rd[] = {0.0, 0.0, 0.0916858};
        //double r=master.eri.ERI(l1, l3, l2, l4, a1, a3, a2, a4, ra, rc, rb, rd);
        double r1 = master.eri.ERI(l1, l3, l2, l4, a1, a3, a2, a4, ra, rc, rb, rd);
    }

    public static void main(String[] args) {
        test_integral run = new test_integral();
        // run.rapat();
        // run.EK();
        // run.NAI();
        Date date_start = new Date();
        SimpleDateFormat ft = new SimpleDateFormat("hh:mm:ss.SSSSSSS");
        String start = ft.format(date_start);
        //for (int i = 0; i < 222111; i++) {
        run.ERI();

        Date date_stop = new Date();
        String stop = ft.format(date_stop);
        System.out.print(start + " " + stop);
    }

}
