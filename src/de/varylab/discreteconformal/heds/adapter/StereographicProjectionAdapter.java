package de.varylab.discreteconformal.heds.adapter;

import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.math.ComplexUtility;

@Position
public class StereographicProjectionAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

	public StereographicProjectionAdapter() {
		super(CoVertex.class, null, null, double[].class, true, false);
	}

	@Override
	public double[] getVertexValue(CoVertex v, AdapterSet a) {
		double[] pos = a.get(CoPositionAdapter.class, Position.class, v, double[].class);
		Complex s = ComplexUtility.stereographic(pos);
		return new double[] {s.re, s.im, 0.0, 1.0};
	}
	
	
	@Override
	public double getPriority() {
		return 100;
	}
	
}
