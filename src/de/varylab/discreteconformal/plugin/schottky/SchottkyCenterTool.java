package de.varylab.discreteconformal.plugin.schottky;

import de.jtem.java2dx.modelling.DraggablePoint2D;

public class SchottkyCenterTool extends DraggablePoint2D {

	private static final long serialVersionUID = 1L;

	public SchottkyCenterTool() {
		setName("Base Point");
		setModelType("schottky base point");
		setText("0");
	}
	
	@Override
	public SchottkyCenterTool clone() {
		return new SchottkyCenterTool();
	}
	
}
