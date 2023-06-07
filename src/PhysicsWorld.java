import java.util.ArrayList;

public class PhysicsWorld {
	public static double gravity = -9.81;
	PhysicsObject box;
	ArrayList<Wall> walls;
	double[] debug = new double[10];
	boolean debugStop;
	
	public PhysicsWorld() {
		box = Shapes.square(1);
		box.setPos(0, 7);
		box.setRot(Math.PI/4 - 0.001);
		box.setVel(0, 0, 0);
		box.setMoI(1.0/6.0);
		box.setCE(0.9);
		walls = new ArrayList<Wall>();
		walls.add(new Wall(new Pair(0, 1), new Pair(0, 2), 0.9));
		walls.add(new Wall(new Pair(1, 0), new Pair(-2, 2), 0.9));
		walls.add(new Wall(new Pair(-1, 0), new Pair(2, 2), 0.9));
		walls.add(new Wall(new Pair(0, -1), new Pair(0, 12), 0.9));
		debugStop = false;
	}
	
	public void step(double dt) {
		box.step(dt, gravity);
		ArrayList<Pair> temp = box.getBounds();
		for (int i = 0; i < temp.size(); i++) {
//			if (temp.get(i).y <= 0) {
//				double r = temp.get(i).x - box.pos.x;
//				
////				Coefficients for quadratic formula:
////				
////				au^2 + bu + c = 0
////				
////				a = 1/(2*mass) + r^2/(2*MoI)
////				b = v.y + r*vrot
////				c = 0
////				
////					Quadratic formula for impulse for elastic collision
////					in the simplified form, u = -b/a:
////				
//				double u = -(box.v.y + r*box.vrot)/(1/(2*box.mass) + r*r/(2*box.MoI));
//				
//				if (u > 0) {
//					double dv = u / box.mass;
//						// 
//					double dw = u * r / box.MoI;
//					box.addV(0, dv, dw);
//				}
//			}
			
//			This is an incomplete implementation of ground collisions.
//			
//			if (temp.get(i).y <= 0) {
////				Pair impulse_vec = impulse * unit vector dir
////				double coeff_friction
////				double change_in_angular_momentum = impulse_vec cross R
////				double change_in_trans_momentum = impulse_vec
////				double dv = impulse_vec / mass
////				double dw = impulse_vec cross R / MoI
////				double change_in_trans_energy = mass * dv^2 / 2 + mass * v dot dv
////				double change_in_rot_energy = MoI * dw^2 / 2 + MoI * w * dw
////				change_in_energy = impulse^2/(2 * mass) + 
//				double r = temp.get(i).x - box.pos.x;
////				double usquareCoeff = 1/(2*m) + r^2/(2*MoI);
////				double uCoeff = v.y + r*vrot;
//				double u = -(box.v.y + r*box.vrot)/(1/(2*box.mass) + r*r/(2*box.MoI));
//					// Quadratic equation for force*time for elastic collision where force
//					// is infinite and time is infinitesimal, with finite product of the two.
//				
//				if (u > 0) {
//					double dv = u / box.mass;
//					double dw = u * r / box.MoI;
//					box.addV(0, dv, dw);
//				}
//			}
			
//			for (int j = 0; j < walls.size(); j++) {
//				if (walls.get(j).collides(temp.get(i))) {
//					double rx = temp.get(i).x - box.pos.x;
//					double ry = temp.get(i).y - box.pos.y;
//					
//					
////					Currently, the frictional coefficient is infinite, meaning there is never any slippage
//					double slipConst = 0;
//					
//					double a = slipConst*walls.get(j).Normal.y - box.v.x + ry*box.vrot;
//					double b = 1/box.mass + ry*ry/box.MoI;
//					double c = -rx*ry/box.MoI;
//					double d = slipConst*walls.get(j).Normal.x + box.v.y + rx*box.vrot;
//					double e = rx*ry/box.MoI;
//					double f = -1/box.mass - rx*rx/box.MoI;
//					
////					Solving the system of equations
////					a = b*impulsex + c*impulsey
////					and
////					d = e*impulsex + f*impulsey
////					
////					The collision is completely plastic.
//					double impulsex = (c*d - a*f)/(c*e - b*f);
//					double impulsey = (a*e - b*d)/(c*e - b*f);
//					
//					if (impulsex*walls.get(j).Normal.x + impulsey*walls.get(j).Normal.y > 0) {
//						double dvx = impulsex / box.mass;
//						double dvy = impulsey / box.mass;
//						double dw = (impulsey*rx - impulsex*ry) / box.MoI;
//						box.addV(dvx, dvy, dw);
//					}
//				}
//			}
			
			for (int j = 0; j < walls.size(); j++) {
				double rx = temp.get(i).x - box.pos.x;
				double ry = temp.get(i).y - box.pos.y;
				double nx = walls.get(j).Normal.x;
				double ny = walls.get(j).Normal.y;
				double rvx = (box.v.x) - (box.vrot)*ry;
				double rvy = (box.v.y) + (box.vrot)*rx;
				if (walls.get(j).collides(temp.get(i)) && rvx*nx + rvy*ny < 0) {
					debug[8] = rvx;
					debug[9] = rvy;
					debugStop = true;
					double ce = box.ce * walls.get(j).ce;
					double m = box.mass;
					double MoI = box.MoI;
					
//					double a = nx;
					double b = -1/m - ry*ry/MoI;
					double c = rx*ry/MoI;
					double d = box.v.x - ry*box.vrot;
//					double e = ny;
//					double f = c;
					double g = -1/m - rx*rx/MoI;
					double h = box.v.y + rx*box.vrot;
					
					double k1 = (1/m + ry*ry/MoI)*nx*nx - 2*nx*ny*rx*ry/MoI + (1/m + rx*rx/MoI)*ny*ny;
					debug[0] = k1;
					double k2 = (1/m + ry*ry/MoI)*nx*ny + (nx*nx - ny*ny)*rx*ry/MoI - (1/m + rx*rx/MoI)*nx*ny;
					double k3 = -box.v.x*nx - box.v.y*ny + nx*ry*box.vrot - ny*rx*box.vrot;
					double k4 = k2/k1;
					double k5 = k3/k1;
					double k6 = (ny*b - nx*c)*nx + (ny*c - nx*g)*ny;
					debug[1] = k6;
					double k7 = (ny*b - nx*c)*ny - (ny*c - nx*g)*nx;
					double k8 = ny*d - nx*h;
					double k9 = k7/k6;
					double k10 = k8/k6;
					
					double rs_const = 2*k10*(box.v.x*nx + box.v.y*ny) + k10*k10*(nx*nx + ny*ny)/m + 2*box.vrot*k10*(rx*ny - ry*nx);
					rs_const += k10*k10*(ny*ny*rx*rx + nx*nx*ry*ry - 2*nx*ny*rx*ry)/MoI;
					
					double ls_const = 2*k5*(box.v.x*nx + box.v.y*ny) + k5*k5*(nx*nx + ny*ny)/m + 2*box.vrot*k5*(rx*ny - ry*nx);
					ls_const += k5*k5*(ny*ny*rx*rx + nx*nx*ry*ry - 2*nx*ny*rx*ry)/MoI;
					ls_const *= (1 - ce);
					
					double rs_coeff = -2*k9*(box.v.x*nx + box.v.y*ny) - 2*k9*k10*(nx*nx + ny*ny)/m - 2*box.vrot*k9*(rx*ny - ry*nx);
					rs_coeff += 2*k10*(-k9*(ny*ny*rx*rx + nx*nx*ry*ry) + nx*ny*(ry*ry - rx*rx))/MoI;
					rs_coeff += 2*k10*(2*k9*nx*ny*rx*ry + rx*ry*(nx*nx - ny*ny))/MoI;
					
					double ls_coeff = -2*k4*(box.v.x*nx + box.v.y*ny) - 2*k4*k5*(nx*nx + ny*ny)/m - 2*box.vrot*k4*(rx*ny - ry*nx);
					ls_coeff += 2*k5*(-k4*(ny*ny*rx*rx + nx*nx*ry*ry) + nx*ny*(ry*ry - rx*rx))/MoI;
					ls_coeff += 2*k5*(2*k4*nx*ny*rx*ry + rx*ry*(nx*nx - ny*ny))/MoI;
					ls_coeff *= (1 - ce);
					
					double rs_sq = k9*k9*(nx*nx + ny*ny)/m + (k9*k9*(ny*ny*rx*rx + nx*nx*ry*ry) + 2*k9*nx*ny*(rx*rx - ry*ry))/MoI;
					rs_sq += rx*ry*(-2*k9*k9*nx*ny + 2*k9*(ny*ny - nx*nx))/MoI;
					
					double ls_sq = k4*k4*(nx*nx + ny*ny)/m + (k4*k4*(ny*ny*rx*rx + nx*nx*ry*ry) + 2*k4*nx*ny*(rx*rx - ry*ry))/MoI;
					ls_sq += rx*ry*(-2*k4*k4*nx*ny + 2*k4*(ny*ny - nx*nx))/MoI;
					ls_sq *= (1 - ce);
					
					double a1 = ls_sq - rs_sq;
					debug[2] = a1;
					double b1 = ls_coeff - rs_coeff;
					double c1 = ls_const - rs_const;
					
					double ipf = (-b1 - Math.sqrt(b1*b1 - 4*a1*c1))/(2*a1);
					
					double ipt = k5 - k4*ipf;
					
					double ipn = k10 - k9*ipf;
					
					double maxE = m*(box.v.x + ipf*ny/m)*(box.v.x + ipf*ny/m);
					maxE += m*(box.v.y - ipf*nx/m)*(box.v.y - ipf*nx/m);
					maxE += MoI*(box.vrot - (ipf*ny*ry + ipf*nx*rx)/MoI)*(box.vrot - (ipf*ny*ry + ipf*nx*rx)/MoI);
					maxE /= 2;
					
					double minE = m*(box.v.x + (ipt*nx + ipf*ny)/m)*(box.v.x + (ipt*nx + ipf*ny)/m);
					minE += m*(box.v.y + (ipt*ny - ipf*nx)/m)*(box.v.y + (ipt*ny - ipf*nx)/m);
					minE += MoI*(box.vrot + ((ipt*ny - ipf*nx)*rx - (ipt*nx + ipf*ny)*ry)/MoI)*(box.vrot + ((ipt*ny - ipf*nx)*rx - (ipt*nx + ipf*ny)*ry)/MoI);
					minE /= 2;
					
					double outE = m*(box.v.x + (ipn*nx + ipf*ny)/m)*(box.v.x + (ipn*nx + ipf*ny)/m);
					outE += m*(box.v.y + (ipn*ny - ipf*nx)/m)*(box.v.y + (ipn*ny - ipf*nx)/m);
					outE += MoI*(box.vrot + ((ipn*ny - ipf*nx)*rx - (ipn*nx + ipf*ny)*ry)/MoI)*(box.vrot + ((ipt*ny - ipf*nx)*rx - (ipt*nx + ipf*ny)*ry)/MoI);
					outE /= 2;
					
					double inE = box.calcEnergy();
					
					debug[7] = 0;
					
					rvx = (box.v.x + (ipn*nx + ipf*ny)/m) - (box.vrot + ((ipn*ny - ipf*nx)*rx - (ipn*nx + ipf*ny)*ry)/MoI)*ry;
					rvy = (box.v.y + (ipn*ny - ipf*nx)/m) + (box.vrot + ((ipn*ny - ipf*nx)*rx - (ipn*nx + ipf*ny)*ry)/MoI)*rx;
					
					if (maxE -1 > inE /*|| outE > maxE*/ || rvx*nx + rvy*ny < 0) {
						ipf = (-b1 + Math.sqrt(b1*b1 - 4*a1*c1))/(2*a1);
						ipn = k10 - k9*ipf;
						debug[7] = 1;
					}
					
					
					double ipx = ipn*nx + ipf*ny;
					double ipy = ipn*ny - ipf*nx;
					
					double dvx = ipx/m;
					double dvy = ipy/m;
					double dw = (ipy*rx - ipx*ry)/MoI;
					
					debug[3] = inE;
					debug[4] = outE;
					debug[5] = minE;
					debug[6] = maxE;
					
					box.addV(dvx, dvy, dw);
				}
			}
		}
	}
}
