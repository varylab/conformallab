package de.varylab.discreteconformal.functional;

import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.halfedgetools.functional.Energy;
import de.jtem.halfedgetools.functional.Functional;
import de.jtem.halfedgetools.functional.Gradient;
import de.jtem.halfedgetools.functional.Hessian;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

public class HyperIdealCirclePatternFunctional implements Functional<CoVertex, CoEdge, CoFace> {

	@Override
	public <HDS extends HalfEdgeDataStructure<CoVertex, CoEdge, CoFace>> void evaluate(
		HDS hds, 
		DomainValue x, 
		Energy E, 
		Gradient G, 
		Hessian H
	) {
		
	}

	@Override
	public boolean hasHessian() {
		return false;
	}

	@Override
	public boolean hasGradient() {
		return true;
	}

	@Override
	public <HDS extends HalfEdgeDataStructure<CoVertex, CoEdge, CoFace>> int getDimension(HDS hds) {
		return hds.numEdges() / 2;
	}

	@Override
	public <HDS extends HalfEdgeDataStructure<CoVertex, CoEdge, CoFace>> int[][] getNonZeroPattern(HDS hds) {
		return null;
	}

}
