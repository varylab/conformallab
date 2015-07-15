package de.varylab.discreteconformal.heds.adapter;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;

import java.util.HashMap;
import java.util.Map;

import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Radius;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

@Radius
public class PseudoConformalStructureAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, Double> {

	private Map<CoEdge, Double>
		lcrPseudo = new HashMap<>();
	
	public PseudoConformalStructureAdapter(Map<CoEdge, Double> lcrPseudo) {
		super(CoVertex.class, null, null, Double.class, true, false);
		this.lcrPseudo = lcrPseudo;
	}
	

	@Override
	public Double getVertexValue(CoVertex v, AdapterSet a) {
		double product = 1.0;
		for (CoEdge e : incomingEdges(v)) {
			Double lcr = lcrPseudo.get(e);
			if (lcr == null) continue;
			product *= lcr;
		}
		return product;
	}
	
	
	@Override
	public double getPriority() {
		return 1;
	}
	
}
