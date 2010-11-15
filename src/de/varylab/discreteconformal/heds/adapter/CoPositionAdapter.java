package de.varylab.discreteconformal.heds.adapter;

import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

@Position
public class CoPositionAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

	public CoPositionAdapter() {
		super(CoVertex.class, null, null, double[].class, true, true);
	}
	
	
	@Override
	public void setVertexValue(CoVertex v, double[] value, AdapterSet a) {
		if (value.length == 3) {
			v.P[0] = value[0];
			v.P[1] = value[1];
			v.P[2] = value[2];
			v.P[3] = 1.0;
		} else {
			v.P = value;
		}
	}
	

	@Override
	public double[] getVertexValue(CoVertex v, AdapterSet a) {
		return v.P;
	}
	
	
	@Override
	public double getPriority() {
		return 1;
	}
	
}
