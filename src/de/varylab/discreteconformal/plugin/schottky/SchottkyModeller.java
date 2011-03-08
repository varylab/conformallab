package de.varylab.discreteconformal.plugin.schottky;

import java.awt.GridLayout;
import java.util.List;

import javax.swing.JFrame;

import de.jtem.java2d.SceneComponent;
import de.jtem.java2d.Viewer2D;
import de.jtem.java2dx.modelling.ModellingTool2D;
import de.jtem.java2dx.modelling.SimpleModeller2D;
import de.jtem.modelling.Modeller;

public class SchottkyModeller extends SimpleModeller2D {

	private SchottkyModelContainer 
		genContainer = new SchottkyModelContainer();
	
	public SchottkyModeller() {
		ModellingTool2D gen = new SchottkyModellingTool();
		gen.setName("schottky generator");
		gen.setModelType("Generator");
		SceneComponent genScene = new SceneComponent();
		Modeller m = getModeller();
		Viewer2D v = getViewer();
		m.addTemplate(gen, genScene);
		m.setModelContainer(genContainer);
		v.setGridEnabled(true);
		v.setTranslateToolEnabled(true);
		setActionTree(m.getViewerMenu());
	}

	
	public List<SchottkyGenerator> getGenerators() {
		return genContainer.getGenerators();
	}
	
	
	public static void main(String[] args) {
		JFrame frame = new JFrame("Schottky Modeller");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setSize(600,600);
		SchottkyModeller modeller = new SchottkyModeller();
		
//		Complex A = new Complex(-1.5, -1.5);
//		Complex B = new Complex(0.0, 0.0);
//		Complex m = new Complex(0.01);
//		A = new Complex(0.23, 0.1);
//		B = new Complex(-0.23, 0.1);
//		m = new Complex(0.005);
//		Complex c = new Complex(0.23, 0.1);
//		SchottkyCircle C1 = new SchottkyCircle(c, 1.0, false);
//		SchottkyGenerator G1 = new SchottkyGenerator(A, B, m, C1);
//		SchottkyModellingTool tool1 = new SchottkyModellingTool(G1); 
//		modeller.getModeller().addTool(tool1, "schottky generator");
//		modeller.getModeller().updateTools();
		
		frame.setLayout(new GridLayout());
		frame.add(modeller.getViewer());
		frame.setVisible(true);
	}

	
}
