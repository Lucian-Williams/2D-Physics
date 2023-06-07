
public class Wall {
	Pair Normal;
	Pair Origin;
	double ce;
	
	public Wall(Pair n, Pair o) {
		Normal = n;
		Origin = o;
		ce = 0.9;
	}
	
	public Wall(Pair n, Pair o, double val) {
		Normal = n;
		Origin = o;
		ce = val;
	}
	
	public boolean collides(Pair p) {
		double offsetx = p.x - Origin.x;
		double offsety = p.y - Origin.y;
		
		return offsetx*Normal.x + offsety *Normal.y <= 0;
	}
}
