package de.varylab.discreteconformal.heds;

public class CustomEdgeInfo {

	public boolean 
		circularHoleEdge = false;
	public double
		phi = Math.PI;
	
	public CustomEdgeInfo() {
	}

	public CustomEdgeInfo(CustomEdgeInfo ei) {
		super();
		this.circularHoleEdge = ei.circularHoleEdge;
		this.phi = ei.phi;
	}
	
}
