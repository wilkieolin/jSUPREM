/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package psuprem;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.distribution.NormalDistribution;
import javax.swing.*;
import org.math.plot.*;
import java.util.ArrayList;

/**
 * Main class to carry out simulations
 * @author Wilkie
 */
public class PSuprem extends Physics{

    /**
     * @param args the command line arguments
     */
    static Node[] points;
    static int NODES;
    static double temp, n_i, e_g;
    public static double Izero_eq, Vzero_eq;
    static Gas ambient;

    public static void main(String[] args) {
        System.out.println("Hello my precious");
        
        //depth is controlled by the number of nodes, with 1 node being 1 nm
        //initialize the grid
        NODES = 2500;
        
        Interface main = new Interface();
        main.show();

        
        gridInit(Math.pow(10,10),Math.pow(10,10));

        set_temp(1273);
        
        /*     ~~~ TEST CODE ~~~
        //preDep(Math.pow(10,20),30*60, Material.PSG);
        preDep(Math.pow(10,20),30*30, Material.BSG);
        //ArrayList<Integer> depths = get_junction_depths();
        //System.out.println(depths.get(0));
        test straggle and range-
        System.out.println(get_dose(Material.BORON));
        System.out.println("RP B: "+get_Rp_B(80.0));
        System.out.println("S B: "+get_S_B(80.0)); 
        System.out.println("RP P: "+get_Rp_P(100.0));
        System.out.println("S P: "+get_S_P(100.0)); 
        
        //graph(plot, xloc, Material.NET_DOPING);
        implant(100, 4*Math.pow(10,15),Material.PHOSPHORUS);
        //graph(Material.BORON);
        //add_oxide(20);
        
        diffuse(3600);
        System.out.println("Mobility : "+points[80].find_mobility()+" cm^2/v*s");
        System.out.println("Mobility : "+points[NODES-1].find_mobility()+" cm^2/v*s");
        System.out.println("Sheet resistance: "+get_Rs()+" ohm/sq");
        System.out.println(get_dose(Material.BORON));
        //graph(Material.NET_DOPING);
        
        overlay_all(); */
        
        System.out.println("Goodbye my sweet");
    }
    public static void gridInit (double B_bg, double P_bg) {
        points = new Node[NODES];    
    
        for (int x=0; x<NODES; x++) { //initialize the grid with background doping 1E10
            points[x] = new Node(B_bg,P_bg,Material.SI);
        }
        set_temp(300);
    }
    
    public static JFrame create_plot () {
        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(600,600);
        frame.setVisible(true);
        
        return frame;
    }
    public static void overlay_all () {
        double[] concs = new double[NODES];
        Plot2DPanel plot = new Plot2DPanel();
        JFrame frame = new JFrame("Boron, phosphorus, and net doping overlay");
        //frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(600,600);
        frame.setContentPane(plot);
        frame.setVisible(true);
        
        double[] xlocs = new double[NODES]; //create array with x-locs for graphing x-axis
        for (int i=0; i<NODES; i++) { xlocs[i]=i; }
        
        
        for (int i=0; i<NODES; i++) { concs[i]=Math.log10(points[i].getB()); } //pull out the concentrations to graph
        plot.addLinePlot("Boron", xlocs, concs);
        
        for (int i=0; i<NODES; i++) { concs[i]=Math.log10(points[i].getP()); } //pull out the concentrations to graph
        plot.addLinePlot("Phosphorus", xlocs, concs);
        
        double conc;
        for (int i=0; i<NODES; i++) { 
            conc = Math.log10(Math.abs(points[i].getNet()));
            if (conc <= 0) {
                conc = 0;
            }
            concs[i] = conc;
        } 
        plot.addLinePlot("Netdoping", xlocs, concs);
    }
    
    public static void graph ( Material newMat) {
        double[] concs = new double[NODES];
        Plot2DPanel plot = new Plot2DPanel();
        JFrame frame = new JFrame(newMat.toString()+" Plot");
        //frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(600,600);
        frame.setContentPane(plot);
        frame.setVisible(true);
        
        double[] xlocs = new double[NODES]; //create array with x-locs for graphing x-axis
        for (int i=0; i<NODES; i++) { xlocs[i]=i; }
        
        switch (newMat) {
            case BORON:
                for (int i=0; i<NODES; i++) { concs[i]=Math.log10(points[i].getB()); } //pull out the concentrations to graph
                plot.addLinePlot("Boron", xlocs, concs);
            break;
            case PHOSPHORUS:
                for (int i=0; i<NODES; i++) { concs[i]=Math.log10(points[i].getP()); } //pull out the concentrations to graph
                plot.addLinePlot("Phosphorus", xlocs, concs);
            break;
            case NET_DOPING:
                double conc;
                for (int i=0; i<NODES; i++) { 
                    conc = Math.log10(Math.abs(points[i].getNet()));
                    if (conc <= 0) {
                        conc = 0;
                    }
                    concs[i] = conc;
                } 
                plot.addLinePlot("Netdoping", xlocs, concs);
            break;
        }
    }
    
    public static void set_temp (double new_temp) {
        temp = new_temp;
        
        Izero_eq = arrhenius(Math.pow(10,27),3.8,temp);
        Vzero_eq = arrhenius(Math.pow(10,23),2.6,temp);
        e_g = E_G0 - ((ALPHA*Math.pow(temp,2))/(temp+BETA));
        
        n_i = 3.1*Math.pow(10,16)*Math.pow(temp,1.5)*Math.exp(-0.603/(BOLTZ_K*temp));
        System.out.println("@"+temp+" K, Izero_eq: "+Izero_eq+" Vzero_eq: "+Vzero_eq+" Bandgap: "+e_g+" Ni: "+n_i);
        
        for (int x=0; x<NODES; x++) {
            points[x].find_D(Izero_eq, Vzero_eq,temp,n_i);
        }
        
    }
    
    public static void add_oxide (int thickness) {
        for (int i = 0; i < thickness; i++ ) {
            points[i].setMat(Material.SIO2);
        }
        System.out.println("Converted "+thickness+" nm to oxide");
    }
    
    public static double get_dose (Material mat) {
        double dose = 0;
        switch (mat) {
            case BORON:
                for (int i = 0; i < NODES; i++) {
                dose += points[i].getB()*Math.pow(10,-7); //sum "slices" of concentration * width
                }
            break;
            case PHOSPHORUS:
                for (int i = 0; i < NODES; i++) {
                dose += points[i].getP()*Math.pow(10,-7); //sum "slices" of concentration * width
                }
            break;
        }
        System.out.println("Dose of "+mat.toString()+": "+dose);
        return dose;       
    }
    
    public static ArrayList<Integer> get_junction_depths () { //returns all junction occurences in an arraylist
        ArrayList<Integer> x_junctions = new ArrayList();
        
        for (int i = 1; i < NODES; i++) { 
            
            if ((points[i-1].getNet() > 0 && points[i].getNet() < 0) || (points[i-1].getNet() < 0 && points[i].getNet() > 0)) { //look for sign changes
                x_junctions.add(i);
                System.out.println("zero found");
            }
        }
        
        return x_junctions;
    }
    
    public static double get_Rs () {
        ArrayList<Integer> x_j = get_junction_depths();
        if (x_j.isEmpty()) { //if there's no junction, calculate Rs for the entire substrate
            x_j.add(NODES-2);
        }
        double sum = 0; 
        double background;
        System.out.println("Finding sheet resistance to "+x_j.get(0)+" nm");
        
        if (points[x_j.get(0)-1].getNet() > 0) { //p-type
            background = Math.abs(points[x_j.get(0)+1].getB()); //find background doping
            for (int i = 0; i < x_j.get(0); i++) {
                sum += (Math.abs(points[i].getB())-background)*(points[i].find_mobility())*Math.pow(10,-7); //find Rs sum term
            }
        } else if (points[x_j.get(0)-1].getNet() < 0) { //n-type
            background = Math.abs(points[x_j.get(0)+1].getP());
            for (int i = 0; i < x_j.get(0); i++) {
                sum += (Math.abs(points[i].getP())-background)*(points[i].find_mobility())*Math.pow(10,-7);
            }
        } else { //intrinsic or counterdoped, calc not supported
            return 0;
        }
        
        System.out.println("Background doping "+background+" cm^-3");
        return 1/(CHARGE*sum);
    }
    
    public static void diffuse (int final_time) {
        double deltaBC, deltaPC, deltaT, dtx, flux;
        deltaT = 1;
        int i;
        
        System.out.println("Diffusing for "+final_time+" seconds");
        for (int time = 0; time<final_time; time++) {
            //take care of the surface, x=0 case using "mirror concentration"
            i = 0;
            
            //BORON
            dtx = (points[i].getBDiff()*deltaT)/(Math.pow(10,-14)); //find  the (D*delta T)/(delta x)^2 term
            if (dtx > 0.5) { dtx = 0.5; }                           //make sure it's <= 0.5
            deltaBC = dtx*(-1*points[i].getB()+points[i+1].getB());    //find the change in C
            points[i].addB(deltaBC);           //set new concentration Cb + delta Cb
            
            //PHOSPHORUS
            dtx = (points[i].getPDiff()*deltaT)/(Math.pow(10,-14));
            if (dtx > 0.5) { dtx = 0.5; } 
            deltaPC = dtx*(-1*points[i].getP()+points[i+1].getP());
            points[i].addP(deltaPC);

            for (i=1; i<NODES-1; i++) { //use the standard formula for the rest of the points
                
                //BORON
                dtx = (points[i].getBDiff()*deltaT)/(Math.pow(10,-14));
                if (dtx > 0.5) { dtx = 0.5; }
                deltaBC = dtx*(points[i-1].getB()-2*points[i].getB()+points[i+1].getB());
                
                if (points[i].getMat()==Material.SIO2 && points[i+1].getMat()==Material.SI) { //check for oxide-si interface **FIX?**
                    flux = V_TRANSFER*(points[i].getB() - points[i+1].getB()/K0_B); //get interface flux
                    System.out.println("Boundary encountered, "+flux+" 1/cm^2*s interface flux");
                    flux = -1*Math.pow(10,-7)*flux/points[i].getBDiff(); //convert flux to change in concentration
                    System.out.println("change in concentration: "+flux);
                    points[i].addB(flux);
                    points[i+1].addB(-1*flux);
                }
                
                points[i].addB(deltaBC);
                
                //PHOSPHORUS
                dtx = (points[i].getPDiff()*deltaT)/(Math.pow(10,-14));
                if (dtx > 0.5) { dtx = 0.5; }
                deltaPC = dtx*(points[i-1].getP()-2*points[i].getP()+points[i+1].getP());
                
                if (points[i].getMat()==Material.SIO2 && points[i+1].getMat()==Material.SI) { //check for oxide-si interface **FIX?**
                    flux = V_TRANSFER*(points[i].getP() - points[i+1].getP()/K0_P); //get interface flux
                    System.out.println("Boundary encountered, "+flux+" 1/cm^2*s interface flux");
                    flux = -1*Math.pow(10,-7)*flux/points[i].getPDiff(); //convert flux to change in concentration
                    System.out.println("change in concentration: "+flux);
                    points[i].addP(flux);
                    points[i+1].addP(-1*flux);
                }
                
                points[i].addP(deltaPC);
            }
            
            i = NODES-1;
            //take care of back surface, also using "mirror concentration"
            
            //BORON
            dtx = (points[i].getBDiff()*deltaT)/(Math.pow(10,-14)); //find  the (D*delta T)/(delta x)^2 term
            if (dtx > 0.5) { dtx = 0.5; }                           //make sure it's <= 0.5
            deltaBC = dtx*(points[i-1].getB()-points[i].getB());    //find the change in C
            points[i].addB(deltaBC);           //set new concentration Cb + delta Cb
            
            //PHOSPHORUS
            dtx = (points[i].getPDiff()*deltaT)/(Math.pow(10,-14));
            if (dtx > 0.5) { dtx = 0.5; }
            deltaPC = dtx*(points[i-1].getP()-points[i].getP());
            points[i].addP(deltaPC);

        }       
       
    }
    
    public static void implant (double energy, double dose, Material mat) {
        NormalDistribution distro;
        double peak_C, straggle, range;
        
        switch (mat) {
            case BORON:
                straggle = get_S_B(energy);
                range = get_Rp_B(energy);
                peak_C = dose/(Math.sqrt(2*Math.PI)*straggle*Math.pow(10,-7));
                for (int i = 0; i < NODES; i++) {
                    points[i].addB(peak_C*Math.exp(-1*Math.pow((i-range),2)/(2*Math.pow(straggle,2))));
                }
            break;
            case PHOSPHORUS:
                straggle = get_S_P(energy);
                range = get_Rp_P(energy);
                peak_C = dose/(Math.sqrt(2*Math.PI)*straggle*Math.pow(10,-7));
                for (int i = 0; i < NODES; i++) {
                    points[i].addP(peak_C*Math.exp(-1*Math.pow((i-range),2)/(2*Math.pow(straggle,2))));
                }
            break;
        }
    }
    
    public static void preDep (double Csurface, double time, Material newMat) {
        switch (newMat) {
            case BSG:
                System.out.println("Boron pre-dep");
                for (int x=0; x<NODES; x++) { //do pre-dep
                    double concB = Csurface*Erf.erfc((x*Math.pow(10,-7))/(2*Math.sqrt(points[x].getBDiff()*time)));
                    //System.out.println(concB+points[x].getB());
                    points[x].addB(concB);
                }
            break;
            case PSG:
                System.out.println("Phosphorus pre-dep");
                for (int x=0; x<NODES; x++) { //do pre-dep
                    double concP = Csurface*Erf.erfc((x*Math.pow(10,-7))/(2*Math.sqrt(points[x].getPDiff()*time)));
                    //System.out.println(concP);
                    points[x].addP(concP);
                }
            break;
        }

    }
    
    public static void gas_flood (Gas new_gas) {
        ambient = new_gas;
    }
    
}
