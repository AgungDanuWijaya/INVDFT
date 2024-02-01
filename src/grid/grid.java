/**
 * Setting grid untuk menghitung Exc dan Vxc
 */
package grid;

import java.util.HashMap;
import java.util.Map;
import function.main_function;
import geo_molecule.data_geo;
import jqc.getdata;

/**
 *
 * @author Richard P. Muller, Agung Danu Wijaya
 */
public class grid {

    public HashMap<Integer, HashMap<Integer, double[]>> points = new HashMap<>();
    grid_atom grd = new grid_atom();

    public grid(String nama, main_function kernel) {
        data_geo geo = kernel.geo.data.get(nama);
        String namaatom[] = geo.atom;
        for (int j = 0; j < namaatom.length; j++) {
            HashMap<Integer, double[]> point = new HashMap<>();
            double batom = geo.numProton[j];
            double[] xyz = geo.R[j];
            point = grd.rungrid(batom, xyz, kernel);
            points.put(j, point);
        }
        becke_reweight_atoms d = new becke_reweight_atoms();
        d.clc_becke_reweight_atoms(kernel, nama, points);
    }

      public double[][] set_gr(main_function kernel, Map<Integer, getdata.datakHF> bfs,
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        int nbf = bfs.size();
        double[][] bfamps = new double[points.get(0).size() * points.size()][nbf];
        for (int j = 0; j < nbf; j++) {
            int index = 0;
            for (int i = 0; i < points.size(); i++) {
                for (int k = 0; k < points.get(i).size(); k++) {
                    double R[] = new double[3];
                    double RW[] = points.get(i).get(k);
                    for (int l = 0; l < R.length; l++) {
                        R[l] = RW[l];
                    }
                    bfamps[index][j] = kernel.gpoint.kinetic(bfs.get(j), R);
                    index++;
                }
            }
        }
        return bfamps;
    }
    public double[][] setbfamps(main_function kernel, Map<Integer, getdata.datakHF> bfs,
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        int nbf = bfs.size();
        double[][] bfamps = new double[points.get(0).size() * points.size()][nbf];
        for (int j = 0; j < nbf; j++) {
            int index = 0;
            for (int i = 0; i < points.size(); i++) {
                for (int k = 0; k < points.get(i).size(); k++) {
                    double R[] = new double[3];
                    double RW[] = points.get(i).get(k);
                    for (int l = 0; l < R.length; l++) {
                        R[l] = RW[l];
                    }
                    bfamps[index][j] = kernel.gpoint.bf(bfs.get(j), R);
                    index++;
                }
            }
        }
        return bfamps;
    }

    public double[][] setbfamps_x(main_function kernel, Map<Integer, getdata.datakHF> bfs,
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        int nbf = bfs.size();
        double[][] bfamps = new double[points.get(0).size() * points.size()][nbf];
        for (int j = 0; j < nbf; j++) {
            int index = 0;
            for (int i = 0; i < points.size(); i++) {
                for (int k = 0; k < points.get(i).size(); k++) {
                    double R[] = new double[3];
                    double RW[] = points.get(i).get(k);
                    for (int l = 0; l < R.length; l++) {
                        R[l] = RW[l];
                    }
                    R[0] += kernel.delta_gama;
                    bfamps[index][j] = kernel.gpoint.bf(bfs.get(j), R);
                    index++;
                }
            }
        }
        return bfamps;
    }

    public double[][] setbfamps_y(main_function kernel, Map<Integer, getdata.datakHF> bfs,
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        int nbf = bfs.size();
        double[][] bfamps = new double[points.get(0).size() * points.size()][nbf];
        for (int j = 0; j < nbf; j++) {
            int index = 0;
            for (int i = 0; i < points.size(); i++) {
                for (int k = 0; k < points.get(i).size(); k++) {
                    double R[] = new double[3];
                    double RW[] = points.get(i).get(k);
                    for (int l = 0; l < R.length; l++) {
                        R[l] = RW[l];
                    }
                    R[1] += kernel.delta_gama;
                    bfamps[index][j] = kernel.gpoint.bf(bfs.get(j), R);
                    index++;
                }
            }
        }
        return bfamps;
    }

    public double[][] setbfamps_z(main_function kernel, Map<Integer, getdata.datakHF> bfs,
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        int nbf = bfs.size();
        double[][] bfamps = new double[points.get(0).size() * points.size()][nbf];
        for (int j = 0; j < nbf; j++) {
            int index = 0;
            for (int i = 0; i < points.size(); i++) {
                for (int k = 0; k < points.get(i).size(); k++) {
                    double R[] = new double[3];
                    double RW[] = points.get(i).get(k);
                    for (int l = 0; l < R.length; l++) {
                        R[l] = RW[l];
                    }
                    R[2] += kernel.delta_gama;
                    bfamps[index][j] = kernel.gpoint.bf(bfs.get(j), R);
                    index++;
                }
            }
        }
        return bfamps;
    }
     public double[][] setbfamps_x_(main_function kernel, Map<Integer, getdata.datakHF> bfs,
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        int nbf = bfs.size();
        double[][] bfamps = new double[points.get(0).size() * points.size()][nbf];
        for (int j = 0; j < nbf; j++) {
            int index = 0;
            for (int i = 0; i < points.size(); i++) {
                for (int k = 0; k < points.get(i).size(); k++) {
                    double R[] = new double[3];
                    double RW[] = points.get(i).get(k);
                    for (int l = 0; l < R.length; l++) {
                        R[l] = RW[l];
                    }
                    R[0] -= kernel.delta_gama;
                    bfamps[index][j] = kernel.gpoint.bf(bfs.get(j), R);
                    index++;
                }
            }
        }
        return bfamps;
    }

    public double[][] setbfamps_y_(main_function kernel, Map<Integer, getdata.datakHF> bfs,
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        int nbf = bfs.size();
        double[][] bfamps = new double[points.get(0).size() * points.size()][nbf];
        for (int j = 0; j < nbf; j++) {
            int index = 0;
            for (int i = 0; i < points.size(); i++) {
                for (int k = 0; k < points.get(i).size(); k++) {
                    double R[] = new double[3];
                    double RW[] = points.get(i).get(k);
                    for (int l = 0; l < R.length; l++) {
                        R[l] = RW[l];
                    }
                    R[1] -= kernel.delta_gama;
                    bfamps[index][j] = kernel.gpoint.bf(bfs.get(j), R);
                    index++;
                }
            }
        }
        return bfamps;
    }

    public double[][] setbfamps_z_(main_function kernel, Map<Integer, getdata.datakHF> bfs,
            HashMap<Integer, HashMap<Integer, double[]>> points) {
        int nbf = bfs.size();
        double[][] bfamps = new double[points.get(0).size() * points.size()][nbf];
        for (int j = 0; j < nbf; j++) {
            int index = 0;
            for (int i = 0; i < points.size(); i++) {
                for (int k = 0; k < points.get(i).size(); k++) {
                    double R[] = new double[3];
                    double RW[] = points.get(i).get(k);
                    for (int l = 0; l < R.length; l++) {
                        R[l] = RW[l];
                    }
                    R[2] -= kernel.delta_gama;
                    bfamps[index][j] = kernel.gpoint.bf(bfs.get(j), R);
                    index++;
                }
            }
        }
        return bfamps;
    }


    public HashMap<Integer, double[]> pointsmap(HashMap<Integer, HashMap<Integer, double[]>> points) {
        HashMap<Integer, double[]> data = new HashMap<>();
        int index = 0;
        for (int i = 0; i < points.size(); i++) {
            for (int k = 0; k < points.get(i).size(); k++) {
                data.put(index, points.get(i).get(k));
                index++;
            }
        }
        return data;
    }

}
