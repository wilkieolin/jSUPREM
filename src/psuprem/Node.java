package psuprem;

/**
 * Class node provides an object which contains the parameters contained in each node of the calculation
 * @author Wilkie
 */

public class Node extends Physics {
    private double boron, phos, net;
    private double Vneg, Vzero, Vplus;
    private double Ineg, Izero, Iplus;
    private double D_boron, D_phos;

    Material mat;
    
    public Node () {
        boron = phos = net = 0;
        setTemp(300);
    }
    public Node (double initBoron, double initPhos, Material initMat) {
        boron = initBoron;
        phos = initPhos;
        net = initBoron - initPhos;
        Izero = Vzero = 0;
        mat = initMat;
    }
    
    public void setB (double conc) {
        if (conc > 0) {
            boron = conc;
            net = phos - boron;
        }
        else System.out.println("bad concentration");
    }
    public void setP (double conc) {
        if (conc > 0) {
            phos = conc;
            net = phos - boron;
        }
        else System.out.println("bad concentration");
    }
    
    public void addB (double conc) { boron += conc; net += conc;}
    public void addP (double conc) { phos += conc; net -= conc;}
    
    public double getB () { return boron; }
    public double getP () { return phos; }
    public double getNet () { return net; }
    public double getVneg () { return Vneg; }
    public double getVzero () { return Vzero; }
    public double getVplus () { return Vplus; }
    public double getIneg () { return Ineg; }
    public double getIzero () { return Izero; }
    public double getIplus () { return Iplus; }
    public double getBDiff () { return D_boron; }
    public double getPDiff () { return D_phos; }
    public Material getMat () { return mat; }
       
    public void setTemp (double temp) {
        //set neutral interstitial and vacancy concentrations
        Izero = arrhenius(Math.pow(10,27),3.8,temp);
        Vzero = arrhenius(Math.pow(10,23),2.6,temp);
    }
    
    public void find_D (double eq_Izero, double eq_Vzero, double temp, double ni) {    
        Izero = eq_Izero;
        Vzero = eq_Vzero;
        
        // n ~ P doping+ ni, p ~ B doping + ni
        // boron = D0 + D+ (p/ni)
        D_boron = arrhenius(B_Di_zero[0],B_Di_zero[1],temp) + arrhenius(B_Di_plus[0],B_Di_plus[1],temp)*((boron+ni)/ni);
        // phos = D0 + D- (n/ni) + D-- (n/ni)^2
        D_phos = arrhenius(P_Di_zero[0],P_Di_zero[1],temp) + arrhenius(P_Di_neg[0],P_Di_neg[1],temp)*((phos+ni)/ni)+arrhenius(P_Di_neg[0],P_Di_2neg[1],temp)*Math.pow(((phos+ni)/ni),2);
        
    }
    
    public double find_mobility () {
        double mobility;
        if (net > 0) { //p-type
            mobility = 54.3 + (406.9)/(1+Math.pow(boron/(2.35*Math.pow(10,17)),0.88));
        } else { //n-type
            mobility = 92 + (1268)/(1+Math.pow(phos/(1.3*Math.pow(10,17)),0.91));
        }
        return mobility;
    }
    
    public void setMat (Material new_mat) {
        mat = new_mat;
    }

}
