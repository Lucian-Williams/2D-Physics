import java.util.ArrayList;

public class PhysicsObject {
	ArrayList<Pair> shape;
	Pair pos; // Position of center of gravity in meters
	Pair v; // Velocity in meters per second
	double rot, vrot; /* Orientation and angular velocity in radians and
						 radians per second respectively */
	double mass; // mass in kilograms
	double MoI; // moment of inertia
	double ce; // coefficient of energy retention on collision
	
	private ArrayList<Pair> bounds;
	
	public PhysicsObject(ArrayList<Pair> srcShape) {
		shape = srcShape;
		bounds = new ArrayList<Pair>();
		for (int i = 0; i < srcShape.size(); i++) {
			bounds.add(new Pair(0, 0));
		}
		pos = new Pair(0, 0);
		rot = 0;
		v = new Pair(0, 0);
		vrot = 0;
		MoI = 1;
		mass = 1;
		ce = 0.9;
	}
	
	public PhysicsObject(ArrayList<Pair> srcShape,
						  double xVal, double yVal, double rotVal) {
		shape = srcShape;
		bounds = new ArrayList<Pair>();
		for (int i = 0; i < srcShape.size(); i++) {
			bounds.add(new Pair(0, 0));
		}
		pos = new Pair(xVal, yVal);
		rot = rotVal;
		v = new Pair(0, 0);
		vrot = 0;
		MoI = 1;
		mass = 1;
		ce = 0.9;
	}
	
	public PhysicsObject(ArrayList<Pair> srcShape,
						  double xVal, double yVal, double rotVal,
						  double vxVal, double vyVal, double vrotVal) {
		shape = srcShape;
		bounds = new ArrayList<Pair>();
		for (int i = 0; i < srcShape.size(); i++) {
			bounds.add(new Pair(0, 0));
		}
		pos = new Pair(xVal, yVal);
		rot = rotVal;
		v = new Pair(vxVal, vyVal);
		vrot = vrotVal;
		MoI = 1;
		mass = 1;
		ce = 0.9;
	}
	
	public void setMoI(double val) {
		MoI = val;
	}
	
	public void setMass(double val) {
		mass = val;
	}
	
	public void setCE(double val) {
		ce = val;
	}
	
	public void setPos(double xVal, double yVal) {
		pos.x = xVal;
		pos.y = yVal;
	}
	
	public void setRot(double rotVal) {
		rot = rotVal;
	}
	
	public void setVel(double vxVal, double vyVal, double vrotVal) {
		v.x = vxVal;
		v.y = vyVal;
		vrot = vrotVal;
	}
	
	public void addV(double vxVal, double vyVal, double vrotVal) {
		v.x += vxVal;
		v.y += vyVal;
		vrot += vrotVal;
	}
	
	public ArrayList<Pair> getBounds() {
		for(int i = 0; i < shape.size(); i++) {
			double x = shape.get(i).x;
			double y = shape.get(i).y;
			bounds.get(i).setPair((x*cos(rot) - y*sin(rot)) + pos.x, (x*sin(rot) + y*cos(rot)) + pos.y);
		}
		return bounds;
	}
	
	public double calcEnergy() {
		return mass*v2()/2 + MoI*vrot*vrot/2;
	}
	
	public double calcEnergy(double ipx, double ipy, double rx, double ry) {
		return mass*((v.x + ipx/mass)*(v.x + ipx/mass) + (v.y + ipy/mass)*(v.y + ipy/mass))/2 + 
				MoI*(vrot + (ipy*rx - ipx*ry)/MoI)*(vrot + (ipy*rx - ipx*ry)/MoI)/2;
	}
	
	public double v() {
		return Math.sqrt(v.x*v.x + v.y*v.y);
	}
	
	public double v2() {
		return v.x*v.x + v.y*v.y;
	}
	
	public static double cos(double val) {
		return Math.cos(val);
	}
	
	public static double sin(double val) {
		return Math.sin(val);
	}
	
	public void step(double dt, double g) {
		pos.x += v.x * dt;
		pos.y += v.y * dt + g * dt * dt/2;
		rot += vrot * dt;
		if (rot >= 0) {
			rot = rot % (2 * Math.PI);
		} else {
			rot = 2*Math.PI - ((-rot) % 2*Math.PI);
		}
		v.y += g * dt;
	}
}

class Shapes {
	public static PhysicsObject square(double s) {
		ArrayList<Pair> temp = new ArrayList<Pair>();
		temp.add(new Pair(-s/2,-s/2));
		temp.add(new Pair(s/2,-s/2));
		temp.add(new Pair(s/2,s/2));
		temp.add(new Pair(-s/2,s/2));
		return new PhysicsObject(temp);
	}
}