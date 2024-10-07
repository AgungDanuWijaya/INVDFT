package solid;

// file: RunShellCommandFromJava.java
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

public class makepw {

    public static void main(String[] args) throws IOException, InterruptedException {

        String command = "make pw";
        String env[] = {};
        //   Process proc = Runtime.getRuntime().exec(command,env,"/home/agung/project/solid/q-e-qe-6.6");
        ProcessBuilder builder = new ProcessBuilder("/bin/bash", "-c", command).directory(new File("/home/agung/project/solid/q-e-qe-6.6"));

        Process proc = builder.start();

        // Read the output
        BufferedReader reader
                = new BufferedReader(new InputStreamReader(proc.getInputStream()));

        String line = "";
        while ((line = reader.readLine()) != null) {
            System.out.print(line + "\n");
        }

        proc.waitFor();

    }
}
