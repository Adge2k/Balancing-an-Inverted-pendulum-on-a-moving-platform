import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display.Arrow;
import org.opensourcephysics.display.Circle;
import org.opensourcephysics.frames.PlotFrame;
import org.opensourcephysics.numerics.ODE;
import org.opensourcephysics.numerics.ODESolver;
import org.opensourcephysics.numerics.RK4;


public class Pendulum extends AbstractSimulation implements ODE , KeyListener{
	PlotFrame frame = new PlotFrame("", "height", "mass spring system (RK4)");
	//PlotFrame trace = new PlotFrame("t","height","history");
	
	Circle mass = new Circle();
	Circle base = new Circle();
	// Change this line to use a different ODE solver
	ODESolver odeSolver =  new RK4(this); //new RK4(this); 
	double xmin, xmax, ymin, ymax;

	// Input variables & parameters
	double g = -9.8; // meters per second square
	Line arm, left, right, top, bottom;
	double m = 10; // kg
	double M;
	double l;
	double I;
	double accelx, accely;
	double platX, platY;
	double theta;
	double state[]= new double[11];  // state of the simulation: height, velocity and time
	
	/**
	 * Constructor for the application.  Sets the properties of the two
	 * display areas.
	 *
	 */
	public Pendulum() {
		/*
		 * frame shows an animation of the bouncing ball
		 */
		//frame.addDrawable(mass);
		//frame.addDrawable(mass2);
		
		/*
		 * trace shows the trajectory of the ball over time
		 */
		/**
		trace.setPreferredMinMax(-1.0,+100.0,-11.0,11.0);
		trace.setConnected(0,true);
		trace.setMarkerShape(0, Dataset.NO_MARKER);*/
	}
	
	/**
	 * do one step of the simulation and display results
	 */
	public void doStep() {
		// step the differential equation
		odeSolver.step();
		// display the current state of the mass
		if(state[8]==Math.PI){
			g = 0;
		}else{
			g = -9.8;
		}
		if(Math.sqrt(sq(state[5])+sq(state[4]))!=l){
			//if(state[8]<Math.PI/2.){
				state[4] = state[0]+(l*cos(state[8]));
				state[5] = (state[1]+(l*sin(state[8])));
			/*}else{
				state[4] = state[0]-(l*cos(state[8]));
				state[5] = (state[1]-(l*sin(state[8])));
			}*/
		}
		arm.setXY(state[0], state[1]);
		arm.setAB(state[4]-state[0], state[5]-state[1]);
		setBox();
		mass.setXY(state[4],state[5]);
		base.setXY(state[0], state[1]);
		
	}
	
	/**
	 * return the state of the simulation
	 */
	public double[] getState() {
		return state ;
	}
	
	/**
	 * compute rates of change for the state variables
	 */
	public void getRate(double state[], double rate[]){
		rate[0] = state[2];
		rate[1] = state[3];
		rate[2] = accelx;
		rate[3] = accely;
		//if(state[8]<Math.PI/2){
		rate[4] = -l*sin(state[8]);		// x-velocity
		rate[5] = l*cos(state[8])*state[9];		// y-velocity
		rate[6] = accelx;		// acceleration of mass in x direction
		rate[7] = 0;		// acceleration of mass in y direction
		rate[8] = state[9];	// Omega
		if(state[8]<Math.PI/2.){
			rate[9] = (g/l)*sin(state[8])+accelx/l;	// Alpha
			rate[4] = l*sin(state[8]);
		}else{
			rate[9] = -(g/l)*sin(state[8])+accelx/l;	// Alpha
		}
		rate[10] = 1;	//time
	}
	
	/**
	 * initialize the state of the simulation
	 */
	public void initialize() {
		// initialize and set the step size for the ODE solver
		frame.clearDrawables();
		odeSolver.initialize(control.getDouble("dt"));
		frame.addKeyListener(this);
		xmin = control.getDouble("x-min");
		xmax = control.getDouble("x-max");
		ymin = control.getDouble("y-min");
		ymax = control.getDouble("y-max");
		frame.setPreferredMinMax(xmin, xmax, ymin, ymax);
		m = control.getDouble("mass");
		l = control.getDouble("l");
		platX = control.getDouble("Platform X");
		platY = control.getDouble("Platform Y");
		I = (4/3)*m*sq(l/2);
		accelx = 0;
		accely = 0;
		theta = Math.PI/2.-0.1;
		state[0] = 0;
		state[1] = 0;
		state[2] = 0;
		state[3] = 0;
		state[4] = state[0]+(l*cos(state[8]));
		state[5] = (state[1]+(l*sin(state[8])));
		state[6] = (sq(state[3])-2*l*state[3]*state[10]*cos(state[8])+sq(l)*sq(state[9]))*cos(state[8]);
		state[7] = (sq(state[3])-2*l*state[3]*state[10]*cos(state[8])+sq(l)*sq(state[9]))*sin(state[8]);
		state[8] = theta;
		state[9] = (state[6]*sin(state[8]))/l;
		state[10] = 0;
		arm = new Line();
		left = new Line();
		right = new Line();
		top = new Line();
		bottom = new Line();
		//setBox();
		frame.addDrawable(arm);
		frame.addDrawable(left);
		frame.addDrawable(right);
		frame.addDrawable(top);
		frame.addDrawable(bottom);
		frame.addDrawable(mass);
		frame.addDrawable(base);
		base.pixRadius = 4;
		arm.setXY(state[0], state[1]);
		arm.setAB(state[4]-state[0], state[5]-state[1]);
		setBox();
		mass.setXY(state[4],state[5]);
		base.setXY(state[0], state[1]);
	}
	
	public void setBox(){
		left.setXY(state[0]-platX/2, state[1]-platY);
		left.setAB(0, platY);
		right.setXY(state[0]+platX/2, state[1]-platY);
		right.setAB(0, platY);
		top.setXY(state[0]-platX/2, state[1]);
		top.setAB(platX, 0);
		bottom.setXY(state[0]-platX/2, state[1]-platY);
		bottom.setAB(platX, 0);
	}
				
	
	/**
	 * sets the simulation parameters to their initial values
	 */
	public void reset() {
		control.setValue("Platform X", 6.0);
		control.setValue("Platform Y", 2.0);
		control.setValue("l",5.0);
		control.setValue("v",0.0);
		control.setValue("dt", 0.01);
		control.setValue("mass", 0.2);
		control.setValue("x-min", -40.0);
		control.setValue("x-max", 40.0);
		control.setValue("y-min", -30.0);
		control.setValue("y-max", 30.0);
		
		initialize();
	}
	
	/**
	 * start the simulation application.
	 * @param args command line parameters (not used)
	 */
	public static void main(String args[]) {
		SimulationControl.createApp(new Pendulum());
	}

	@Override
	public void keyPressed(KeyEvent e) {
		if(e.getKeyCode()==KeyEvent.VK_LEFT){
			state[2] = -3;
			accelx = -15;
			
		// TODO Auto-generated method stub
		}
		
		if(e.getKeyCode()==KeyEvent.VK_RIGHT){
			state[2] = 3;
			accelx = 15;
			
		// TODO Auto-generated method stub
		}
		
	}

	@Override
	public void keyReleased(KeyEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void keyTyped(KeyEvent arg0) {
		// TODO Auto-generated method stub
		
	}
	double sq(double a){
		return Math.pow(a, 2.0);
	}
	double cos(double theta){
		return Math.cos(theta);
	}
	
	double sin(double theta){
		return Math.sin(theta);
	}
}
