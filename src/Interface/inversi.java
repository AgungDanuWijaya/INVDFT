/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package Interface;

import java.io.IOException;
import run.inversi_driver;

/**
 *
 * @author agung
 */
public class inversi {
    public static void main(String[] args) throws ClassNotFoundException, IOException {
        inversi_driver a = new inversi_driver();
        Interface gui=new Interface();
        a.manis(gui);
    }
}
