package de.varylab.discreteconformal.heds;

public class CustomEdgeInfo {

	public boolean 
		holeEdge = false;
	
	public CustomEdgeInfo() {
	}

	public CustomEdgeInfo(CustomEdgeInfo ei) {
		super();
		this.holeEdge = ei.holeEdge;
	}
	
}
