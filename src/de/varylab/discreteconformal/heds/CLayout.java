package de.varylab.discreteconformal.heds;

import java.util.LinkedList;

import no.uib.cipr.matrix.Vector;

public class CLayout {

	
	/**
	 * Do flat layout for a HDS and a scale vector u
	 * @param hds
	 * @param u
	 */
	public static void doLayout(CHDS hds, Vector u) {
		LinkedList<CEdge> edgeList = new LinkedList<CEdge>();
		LinkedList<Double> angleList = new LinkedList<Double>();
		
		edgeList.offer(hds.getEdge(0));
		angleList.offer(0.0);
		
		
		
	}
	
	
}
