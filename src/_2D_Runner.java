import java.awt.Frame;
import java.awt.Toolkit;

import javax.swing.JFrame;


public class _2D_Runner{
	
	/**
	 * 
	 */
	static PhysicsPanel mainPanel;
	
	public static void main(String[] args){
		mainPanel = new PhysicsPanel(Toolkit.getDefaultToolkit().getScreenSize());
		JFrame frame = new JFrame("2D");
		frame.setExtendedState(Frame.MAXIMIZED_BOTH);
		frame.setUndecorated(true);
		frame.addKeyListener(mainPanel);
		//frame.setSize(800, 680);
		frame.setFocusable(true);
		frame.setFocusTraversalKeysEnabled(false);
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(mainPanel);
	}
}
