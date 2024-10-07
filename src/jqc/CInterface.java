/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jqc;

       
import com.sun.jna.Library; 
       
public interface CInterface extends Library 
{  
    public void lda_x(double[]rho,double []exc,double []dexc, int n,int id);
    public void gga_x(double[]rho,double[]sigma,double []exc,double []dexc,double []derho, int n,int id);
      void besel_all (double r[],double out[],int pan,double n,double q);
      public void sum(double[]n1,double []n3, int n2);
      public void mdot(double a[], double b,int i,int j,int k ) ;
      public void time_complex(double a[], double b[], double c[],int a_i,int a_j,int a_k,int b_i,int b_j,int b_k,int c_i,int c_j,int c_k);
      public double simpson(int msh, double func[], double rab[]);
}     