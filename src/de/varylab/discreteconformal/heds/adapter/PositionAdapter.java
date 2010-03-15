package de.varylab.discreteconformal.heds.adapter;

import geom3d.Point;
import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

@Position
public class PositionAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

	public PositionAdapter() {
		super(CoVertex.class, null, null, double[].class, true, true);
	}
	
	
	@Override
	public void setVertexValue(CoVertex v, double[] value, AdapterSet a) {
		v.getPosition().set(value);
	}
	

	@Override
	public double[] getVertexValue(CoVertex v, AdapterSet a) {
		Point p = v.getPosition();
		return new double[] {p.x(), p.y(), p.z(), 1.0};
	}
	
	
	@Override
	public double getPriority() {
		return 1;
	}
	
}
