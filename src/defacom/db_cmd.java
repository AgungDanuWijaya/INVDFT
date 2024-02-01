package defacom;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import function.main_function;
import javax.swing.JOptionPane;

public class db_cmd {

    public String inputdata[][];
    public double output[][];
    public double inputdata_T[];
    public double output_T[];
    public double R[][];
    public int numE[];
    public String long_name;
    public String id;

    public void setdata(main_function kernel, String[] strings, HashMap<String, Integer> hmap) throws ParseException {
        String data = strings[3];
        this.id = strings[0];
        this.long_name = strings[1];
        Object obj = new JSONParser().parse(data);
        JSONObject jsonObject = (JSONObject) obj;
        kernel.spinup = Integer.parseInt(jsonObject.get("Spin_up").toString());
        kernel.spindn = Integer.parseInt(jsonObject.get("Spin_dn").toString());
        JSONArray atom = (JSONArray) jsonObject.get("atom");
        int numE[] = new int[atom.size()];
        Iterator<String> iterator = atom.iterator();
        int yer = 0;
        while (iterator.hasNext()) {
            String name = (String) iterator.next();
            numE[yer++] = hmap.get(name);
        }
        this.numE = numE;
        JSONObject xyz = (JSONObject) jsonObject.get("xyz");
        double R[][] = new double[xyz.size()][];
        for (int i = 1; i <= xyz.size(); i++) {
            JSONArray xyz_c = (JSONArray) xyz.get("" + i + "");
            Iterator<Double> iterator_xyz = xyz_c.iterator();
            double R_i[] = new double[3];
            String temp = "";
            while (iterator_xyz.hasNext()) {
                temp += iterator_xyz.next() + "#";
            }
            String temps[] = temp.split("#");
            for (int j = 0; j < temps.length; j++) {
                R_i[j] = Double.parseDouble(temps[j]);
            }
            R[i - 1] = R_i;
        }
        this.R = R;
    }

    public void get_mole(String name_mole, main_function kernel) throws ParseException {
        HashMap<String, Integer> hmap = kernel.atom_n.name_atom();
        //String query = "SELECT * FROM Quantum.xyz_molecule where name='" + name_mole + "';";
        String query = "SELECT * FROM Quantum.xyz_molecule where Simbol='" + name_mole + "';";
       
       // System.err.println(query);
        if (kernel.a.data_src.equals("MySQL")) {
            dbprocess db = new dbprocess(database.local, kernel);
            try {
                ArrayList<String[]> list = db.getArrayData(query);
                for (String[] strings : list) {
                    setdata(kernel, strings, hmap);
                }
                if (list.size() == 0) {
                    JOptionPane.showMessageDialog(null, "Data tidak ditemukan di database !");
                }
            } catch (Exception e) {

            }
        } else {
            String dat[] = kernel.a.data_geo.split(";");
            String strings[] = {"1", dat[0], dat[1], dat[2]};
            setdata(kernel, strings, hmap);
        }
    }
}
