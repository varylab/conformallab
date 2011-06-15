package de.varylab.discreteconformal.heds;

public class CustomEdgeInfo {

	public boolean 
		circularHoleEdge = false;
	
	public CustomEdgeInfo() {
	}

	public CustomEdgeInfo(CustomEdgeInfo ei) {
		super();
		this.circularHoleEdge = ei.circularHoleEdge;
	}
	
}
