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
		super(CoVertex.class, null, CoFace.class, double[].class, true, true);
	}
	
	@Override
	public void setVertexValue(CoVertex v, double[] value, AdapterSet a) {
		if (value.length == 2) {
			v.P[0] = value[0];
			v.P[1] = value[1];
			v.P[2] = 0.0;
			v.P[3] = 1.0;
		} else 
		if (value.length == 3) {
			v.P[0] = value[0];
			v.P[1] = value[1];
			v.P[2] = value[2];
			v.P[3] = 1.0;
		} else 
		if (value.length == 4) {
			System.arraycopy(value, 0, v.P, 0, 4);
		} else {
			throw new IllegalArgumentException("invalid dimension in set vertex value of CoVertex");
		}
	}
	
	@Override
	public void setFaceValue(CoFace f, double[] value, AdapterSet a) {
		if (value.length == 2) {
			f.P[0] = value[0];
			f.P[1] = value[1];
			f.P[2] = 0.0;
			f.P[3] = 1.0;
		} else 
		if (value.length == 3) {
			f.P[0] = value[0];
			f.P[1] = value[1];
			f.P[2] = value[2];
			f.P[3] = 1.0;
		} else
		if (value.length == 4) {
			System.arraycopy(value, 0, f.P, 0, 4);
		} else {
			throw new IllegalArgumentException("invalid dimension in set vertex value of CoFace");
		}
	}
	

	@Override
	public double[] getVertexValue(CoVertex v, AdapterSet a) {
		return v.P;
	}
	
	@Override
	public double[] getFaceValue(CoFace f, AdapterSet a) {
		return f.P;
	}
	
	@Override
	public double getPriority() {
		return 1;
	}
	
}
