package de.varylab.discreteconformal.plugin.schottky;

import java.awt.GridLayout;
import java.util.List;

import javax.swing.JFrame;

import de.jtem.java2d.SceneComponent;
import de.jtem.java2d.Viewer2D;
import de.jtem.java2dx.modelling.ModellingTool2D;
import de.jtem.java2dx.modelling.SimpleModeller2D;
import de.jtem.mfc.field.Complex;
import de.jtem.modelling.ActionTree;
import de.jtem.modelling.Modeller;

public class SchottkyModeller extends SimpleModeller2D {

	private SchottkyModelContainer 
		modelContainer = new SchottkyModelContainer();
	
	public SchottkyModeller() {
		Modeller m = getModeller();
		Viewer2D v = getViewer();
		
		// generator template
		ModellingTool2D generatorTemplate = new SchottkyGeneratorTool();
		SceneComponent genScene = new SceneComponent();
		m.addTemplate(generatorTemplate, genScene);
		
		m.setModelContainer(modelContainer);
		v.setGridEnabled(true);
		v.setTranslateToolEnabled(true);
		v.setScaleToolEnabled(true);
		ActionTree menu = m.getViewerMenu();
//		menu.add(m.getFileMenu());
//		menu.add(m.getEditMenu());
//		menu.add(m.getInspectMenu());
//		setActionTree(menu);
		setActionTree(menu);
		
		m.addTool(new SchottkyCenterTool(), null);
		m.addTool(new SchottkyGeneratorTool(), "Generator");
		m.updateTools();
	}

	
	public List<SchottkyGenerator> getGenerators() {
		return modelContainer.getGenerators();
	}
	
	public void setGenerators(List<SchottkyGenerator> generators) {
		modelContainer.setGenerators(generators);
		getModeller().updateTools();
		getViewer().repaint();
	}
	
	public void removeGenerator(SchottkyGenerator g) {
		modelContainer.removeModel(g, "schottky generator");
	}
	
	public Complex getBasePoint() {
		return modelContainer.getBasePoint();
	}
	
	public SchottkyModelContainer getModelContainer() {
		return modelContainer;
	}
	
	
	public static void main(String[] args) {
		JFrame frame = new JFrame("Schottky Modeller");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setSize(600,600);
		SchottkyModeller modeller = new SchottkyModeller();
		frame.setLayout(new GridLayout());
		frame.add(modeller.getViewer());
		frame.setVisible(true);
	}

	
}
