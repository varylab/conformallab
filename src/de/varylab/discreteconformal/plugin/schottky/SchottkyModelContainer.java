package de.varylab.discreteconformal.plugin.schottky;

import java.util.LinkedList;
import java.util.List;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import de.jtem.modelling.ModelContainer;

public class SchottkyModelContainer implements ModelContainer {

	private List<SchottkyGenerator>
		generators = new LinkedList<SchottkyGenerator>();
	private LinkedList<ChangeListener>
		listeners = new LinkedList<ChangeListener>();

	public List<SchottkyGenerator> getGenerators() {
		return generators;
	}
	
	@Override
	public boolean addModel(Object model, String modelType) {
		generators.add((SchottkyGenerator)model);
		fireStateChanged(this);
		return true;
	}

	@Override
	public void removeModel(Object model, String modelType) {
		generators.remove(model);
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
