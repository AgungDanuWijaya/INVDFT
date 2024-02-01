/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package Interface;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import run.grad_driver;
import run.inversi_driver;

/**
 *
 * @author agung
 */
public class cluster {

    public static void main(String[] args) throws ClassNotFoundException, IOException {
        grad_driver a = new grad_driver();
        Interface gui = new Interface();
        gui.tugas="1";
        a.manis(gui);

    }
}
