package geo_molecule;

public class data_geo {

    public String[] atom; // atom's name
    public int[] numProton; // number electron or proton each atom
    public double[][] R; // atom position

    public data_geo(int[] numElect, String[] atom, double[][] R) {
        this.numProton = numElect;
        this.atom = atom;
        this.R = R;
    }
}
