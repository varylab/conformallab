package de.varylab.discreteconformal.heds;

import de.jtem.halfedge.HalfEdgeDataStructure;

public class CHDS extends HalfEdgeDataStructure<CVertex, CEdge, CFace> {

	public CHDS() {
		super(CVertex.class, CEdge.class, CFace.class);
	}

	
	public void prepareData(double[] theta) {
		for (CEdge e : getEdges()) {
			e.prepareData();
		}
		for (CVertex v : getVertices()) {
			v.setTheta(theta[v.getIndex()]);
		}
	}
	
	
	public void conformalEnergy(double[] u, double[] E, double[] grad, double[][] hess) {
		for (CVertex v : getVertices())
			v.setU(u[v.getIndex()]);
		// Energy
		E[0] = 0.0;
		for (CFace t : getFaces())
			E[0] += t.getFaceEneergy();
		for (CVertex v : getVertices())
			E[0] += v.getVertexEnergy();
		// Gradient
		
	}
	
	
	
	
	
}
