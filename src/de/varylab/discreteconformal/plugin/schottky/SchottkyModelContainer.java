package de.varylab.discreteconformal.plugin.schottky;

import java.util.LinkedList;
import java.util.List;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import de.jtem.java2dx.Point2DDouble;
import de.jtem.mfc.field.Complex;
import de.jtem.modelling.ModelContainer;

public class SchottkyModelContainer implements ModelContainer {

	private List<SchottkyGenerator>
		generators = new LinkedList<SchottkyGenerator>();
	private Point2DDouble
		basePoint = new Point2DDouble();
	private LinkedList<ChangeListener>
		listeners = new LinkedList<ChangeListener>();

	public List<SchottkyGenerator> getGenerators() {
		return generators;
	}
	public void setGenerators(List<SchottkyGenerator> generators) {
		this.generators = generators;
		fireStateChanged(this);
	}
	
	public Complex getBasePoint() {
		return new Complex(basePoint.x, basePoint.y);
	}
	public void setBasePoint(Point2DDouble basePoint) {
		this.basePoint = basePoint;
		fireStateChanged(this);
	}
	
	@Override
	public boolean addModel(Object model, String modelType) {
		if (modelType == null) return false;
		if (modelType.equals("schottky generator")) {
			generators.add((SchottkyGenerator)model);
		}
		if (modelType.equals("schottky base point")) {
			basePoint = (Point2DDouble)model;
		}
		fireStateChanged(this);
		return true;
	}

	@Override
	public void removeModel(Object model, String modelType) {
		if (modelType == null) return;
		if (modelType.equals("schottky generator")) {
			generators.remove(model);
		}
		if (modelType.equals("schottky base point")) {
//			throw new RuntimeException("base point cannot be removed");
		}
		fireStateChanged(this);
	}

	public void fireStateChanged(Object o) {
		for (ChangeListener l : listeners) {
			l.stateChanged(new ChangeEvent(o));
		}
	}
	
	@Override
	public boolean modelChanged(Object model, String modelType) {
		fireStateChanged(this);
		return true;
	}

	@Override
	public void addChangeListener(ChangeListener l) {
		listeners.add(l);
	}

	@Override
	public void removeChangeListener(ChangeListener l) {
		listeners.remove(l);
	}

}
