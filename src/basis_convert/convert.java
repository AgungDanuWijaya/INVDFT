/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package basis_convert;

/**
 *
 * @author agung
 */
public class convert {

    
    public static void main(String[] args) {
     
    String hj=new data().basis;
    String da[]=hj.split("#");
    //    System.out.println(da[1]);
       int index=da[1].indexOf("S");
        int index_=index;
       while (index >= 0) {
   
    index = da[1].indexOf("S", index + 1);
           System.out.println(index);
   //  System.out.println(da[1].substring(index_, index));
     index_=index;
}
    }
    
}
