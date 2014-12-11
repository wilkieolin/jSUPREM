/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package psuprem;

/**
 * Provides class with useful constants and functions
 * @author Wilkie
 */
public class Physics {
    
    final static double BOLTZ_K = 8.617332478*Math.pow(10,-5);
    final static double N_C = 2.8*Math.pow(10,19);
    final static double N_V = 1.04*Math.pow(10,19);
    final static double E_G0 = 1.17;
    final static double BETA = 636;
    final static double ALPHA = 4.73*Math.pow(10,-4);
    final static double CHARGE = 1.602*Math.pow(10,-19);
    
    final static double[] P_Di_zero = {3.85, 3.66};
    final static double[] P_Di_plus = {0, 0};
    final static double[] P_Di_neg = {4.44, 4.0};
    final static double[] P_Di_2neg = {44.2, 4.37};
    
    final static double[] B_Di_zero = {0.05, 3.5};
    final static double[] B_Di_plus = {0.95, 3.5};
    final static double[] B_Di_neg = {0, 0};
    final static double[] B_Di_2neg = {0, 0};   
    
    final static double K0_B = 0.3;
    final static double K0_P = 30;
    final static double V_TRANSFER = 1.55*Math.pow(10,-7);
    
    //final double E_G = 1.12;
    
    public static double arrhenius (double k0, double ea, double temp) {
        return (k0)*Math.exp((-1*ea)/(BOLTZ_K*temp));
    }
    
    public static double get_Rp_B (double energy) {
        return -3.60*Math.pow(10,-3)*Math.pow(energy,2)+3.22*energy+2.4;
    }
    
    public static double get_Rp_P (double energy) {
        return 2.95*Math.pow(10,-4)*Math.pow(energy,2)+1.192*energy+0.5266;
    }
    
    public static double get_S_B (double energy) {
        return -2.3*Math.pow(10,-3)*Math.pow(energy,2)+0.8453*energy+7.126;
    }
    
    public static double get_S_P (double energy) {
        return -6.2*Math.pow(10,-4)*Math.pow(energy,2)+0.498*energy+0.898;
    }
    
}
