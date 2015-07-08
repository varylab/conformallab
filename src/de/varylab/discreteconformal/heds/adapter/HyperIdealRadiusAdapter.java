package de.varylab.discreteconformal.heds.adapter;

import static de.varylab.discreteconformal.math.MathUtility.arsinh;
import static java.lang.Math.sinh;
import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Radius;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

@Radius
public class HyperIdealRadiusAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, Double> {

	private double[] 
		solution = null;
	
	public HyperIdealRadiusAdapter(double[] solution) {
		super(CoVertex.class, null, CoFace.class, Double.class, true, false);
		this.solution = solution;
	}
	
	@Override
	public Double getVertexValue(CoVertex v, AdapterSet a) {
		int index = v.getSolverIndex();
		if (index < 0) {
			return 0.0;
		} else {
			double b = solution[v.getSolverIndex()];
			return arsinh(1.0 / sinh(b));
		}
	}
	
	@Override
	public Double getFaceValue(CoFace f, AdapterSet a) {
		return 1.0;
	}
	
	@Override
	public double getPriority() {
		return 1;
	}
	
}
