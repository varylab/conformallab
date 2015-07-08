package de.varylab.discreteconformal.plugin;

import static java.awt.RenderingHints.KEY_ANTIALIASING;
import static java.awt.RenderingHints.VALUE_ANTIALIAS_ON;

import java.awt.BasicStroke;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Path2D;

import javax.swing.JFrame;

public class VMCrashReproduction extends JFrame {

	private static final long serialVersionUID = 1L;

	public VMCrashReproduction() {
		setSize(1000, 1000);
		setDefaultCloseOperation(EXIT_ON_CLOSE);
	}
	
	@Override
	public void paint(Graphics g) {
		super.paint(g);
		Graphics2D g2d = (Graphics2D)g;
		Path2D path = new Path2D.Float();
		path.append(new Ellipse2D.Double(0.2014256547853074, 0.4506075538048867, 0.2779263027163357, 0.2779263027163357), false);
		path.append(new Ellipse2D.Double(0.5955270094527079, 0.3041804257465449, 0.18761279304214137, 0.18761279304214137), false);
		path.append(new Ellipse2D.Double(-0.26794919129277944, -0.2679491912927795, 0.5358983825855589, 0.5358983825855589), false);
		path.append(new Ellipse2D.Double(-0.4793519570516376, 0.45060755990094614, 0.2779262965084575, 0.2779262965084575), false);
		
		g2d.setRenderingHint(KEY_ANTIALIASING, VALUE_ANTIALIAS_ON);
		g2d.setStroke(new BasicStroke(.01f));
		g2d.translate(100, 100);
		g2d.scale(100.0, 100.0);
		g2d.draw(path);
	}
	
	public static void main(String[] args) {
		new VMCrashReproduction().setVisible(true);
	}
	
}
