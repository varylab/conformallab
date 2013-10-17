package de.varylab.discreteconformal.plugin.schottky;

import java.awt.GridLayout;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JFrame;

import de.jtem.java2d.SceneComponent;
import de.jtem.java2dx.modelling.ModellingTool2D;
import de.jtem.java2dx.modelling.SimpleModeller2D;
import de.jtem.mfc.field.Complex;
import de.jtem.modelling.Modeller;
import de.jtem.modelling.ModellingTool;

public class SchottkyModeller extends SimpleModeller2D {

	private static final String
		TEMPLATE_GENERATOR = "Generator";
	
	public SchottkyModeller() {
		Modeller m = getModeller();
		getModeller().setModelContainer(new SchottkyModelContainer());
		m.clear();
		ModellingTool2D generatorTemplate = new SchottkyGeneratorTool();
		SceneComponent generatorScene = new SceneComponent();
		m.addTemplate(generatorTemplate, generatorScene);
		reset();
		createDefaultData();
	}
	
	public List<SchottkyGenerator> getGenerators() {
		return getModelContainer().getGenerators();
	}
	
	public void setGenerators(List<SchottkyGenerator> generators) {
		reset();
		for (SchottkyGenerator s : generators) {
			SchottkyGeneratorTool tool = new SchottkyGeneratorTool(s);
			getModeller().addTool(tool, TEMPLATE_GENERATOR);
		}
	}
	
	public void reset() {
		getModeller().clear();
	}
	
	public void createDefaultData() {
		Modeller m = getModeller();
		m.addTool(new SchottkyCenterTool(), null);
		m.addTool(new SchottkyGeneratorTool(), TEMPLATE_GENERATOR);
		m.updateTools();
		setActionTree(m.getViewerMenu());
	}
	
	private List<SchottkyGeneratorTool> getGeneratorTools() {
		List<SchottkyGeneratorTool> tools = new LinkedList<SchottkyGeneratorTool>();
		for (ModellingTool tool : getModeller().getToolToTemplate().keySet()) {
			if (tool instanceof SchottkyGeneratorTool) {
				SchottkyGeneratorTool st = (SchottkyGeneratorTool)tool;
				tools.add(st);
			}
		}
		return tools;
	}
	
	
	public void removeGenerator(SchottkyGenerator g) {
		for (SchottkyGeneratorTool tool : getGeneratorTools()) {
			if (tool.getModel() == g) {
				getModeller().removeTool(tool);
			}
		}
	}
	
	public Complex getBasePoint() {
		return getModelContainer().getBasePoint();
	}
	
	public SchottkyModelContainer getModelContainer() {
		return (SchottkyModelContainer)getModeller().getModelContainer();
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
