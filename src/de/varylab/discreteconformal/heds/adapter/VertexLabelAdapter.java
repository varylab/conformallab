package de.varylab.discreteconformal.heds.adapter;

import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Label;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

@Label
public class VertexLabelAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, String> {
	
	public VertexLabelAdapter() {
		super(CoVertex.class, null, null, String.class, true, false);
	}
	
	@Override
	public String getVertexValue(CoVertex v, AdapterSet a) {
		return v.getIndex() + "";
	}
	
	@Override
	public double getPriority() {
		return 1;
	}
	
}
