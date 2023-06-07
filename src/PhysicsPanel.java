import java.util.ArrayList;
import java.awt.Graphics;
import java.awt.Polygon;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

import javax.swing.Timer;
import javax.swing.JPanel;

public class PhysicsPanel extends JPanel implements ActionListener, KeyListener{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private PhysicsWorld world;
	private Dimension dim;
	private double scale;
	private Pair camera;
	private boolean paused;

	Timer tm = new Timer(1, this);

	public PhysicsPanel(Dimension panelSize) {
		paused = true;
		world = new PhysicsWorld();
		dim = panelSize;
		scale = 50;
		camera = new Pair(0,0);
		tm.start();
	}
	
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		
		camera.x = world.box.pos.x;
		camera.y = world.box.pos.y;
		
		int cx = scaleX(camera.x);
		int cy = scaleY(camera.y);
		
		Polygon p = new Polygon();
		ArrayList<Pair> temp = world.box.getBounds();
		
		for (int i = 0; i < temp.size(); i++) {
			p.addPoint(mapX(temp.get(i).x) - cx, mapY(temp.get(i).y) - cy);
		}
		
		g.fillPolygon(p);
		
		for (int i = 0; i < world.walls.size(); i++) {
			
			int x = mapX(world.walls.get(i).Origin.x);
			int y = mapY(world.walls.get(i).Origin.y);
			
			p = new Polygon();
			
			if (world.walls.get(i).Normal.y != 0) {
				int rightHeight = (int)Math.round(world.walls.get(i).Normal.x*(dim.width - x)/world.walls.get(i).Normal.y) + y;
				int leftHeight = (int)Math.round(-world.walls.get(i).Normal.x*x/world.walls.get(i).Normal.y) + y;
				
				if (world.walls.get(i).Normal.y > 0) {
					p.addPoint(dim.width - cx, rightHeight - cy);
					p.addPoint(dim.width - cx, dim.height - cy);
					p.addPoint(0 - cx, dim.height - cy);
					p.addPoint(0 - cx, leftHeight - cy);
				} else {
					p.addPoint(dim.width - cx, rightHeight - cy);
					p.addPoint(0 - cx, leftHeight - cy);
					p.addPoint(0 - cx, 0 - cy);
					p.addPoint(dim.width - cx, 0 - cy);
				}
			} else if (world.walls.get(i).Normal.x > 0) {
				p.addPoint(x - cx, 0 - cy);
				p.addPoint(x - cx, dim.height - cy);
				p.addPoint(0 - cx, dim.height - cy);
				p.addPoint(0 - cx, 0 - cy);
			} else {
				p.addPoint(x - cx, 0 - cy);
				p.addPoint(dim.width - cx, 0 - cy);
				p.addPoint(dim.width - cx, dim.height - cy);
				p.addPoint(x - cx, dim.height - cy);
			}
			
			g.fillPolygon(p);
		}
		
		
//		g.fillRect(0, dim.height - 100, dim.width, dim.height);
		
		g.setColor(Color.RED);
		double kE = world.box.calcEnergy();
		double pE = -world.box.mass * PhysicsWorld.gravity * world.box.pos.y;
		g.drawString("" + (kE + pE), 100, 100);
		
		for (int i = 0; i < world.debug.length; i++) {
			g.drawString("" + world.debug[i], 100, 150 + i*15);
		}
	}
	
	public void actionPerformed(ActionEvent arg0) {
		if (!paused && !world.debugStop)
			for (int i = 0; i < 100; i++) {
				world.step(0.0001);
				if (world.debugStop)
					i += 100;
			}
		
		repaint();
	}
	
	private int scaleX(double val) {
		return (int)Math.round(scale*val);
	}
	
	private int scaleY(double val) {
		return -(int)Math.round(scale*val);
	}
	
	private int mapX(double val) {
		return (int)Math.round(scale*val) + dim.width/2;
	}
	
	private int mapY(double val) {
		return (int)Math.round(-scale*val) + dim.height/2;
	}

	public void keyPressed(KeyEvent e) {
		
		int c = e.getKeyCode();
		
		if (c == KeyEvent.VK_SPACE) {
			if (world.debugStop)
				world.debugStop = false;
			else
				paused = !paused;
		}
		
		if (c == KeyEvent.VK_ESCAPE) {
			System.exit(0);
		}
	}

	public void keyReleased(KeyEvent e) {
		
	}

	public void keyTyped(KeyEvent e) {
		
	}
}
