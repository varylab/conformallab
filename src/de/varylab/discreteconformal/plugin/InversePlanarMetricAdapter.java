package de.varylab.discreteconformal.plugin;

import static java.lang.Math.exp;
import no.uib.cipr.matrix.Vector;
import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Node;
import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.generic.Position3d;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

@Length
public class InversePlanarMetricAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, Double> {

	private Vector
		u = null;
	
	public InversePlanarMetricAdapter(Vector u) {
		super(CoVertex.class, CoEdge.class, CoFace.class, Double.class, true, false);
		this.u = u;
	}
	
	/**
	 * Calculate the edge inverse length for the flat metric
	 * @param e
	 * @param u
	 * @return the new edge length
	 */
	public Double getNewLength(CoEdge e, double length, Vector u) {
		CoVertex v1 = e.getStartVertex();
		CoVertex v2 = e.getTargetVertex();
		Double u1 = v1.getSolverIndex() >= 0 ? u.get(v1.getSolverIndex()) : 0.0; 
		Double u2 = v2.getSolverIndex() >= 0 ? u.get(v2.getSolverIndex()) : 0.0;
		return exp(Math.log(length) - u1 - u2);
	}
	
	
	@Override
	public Double getEdgeValue(CoEdge e, AdapterSet a) {
		double[] s = a.getD(Position3d.class, e.getStartVertex());
		double[] t = a.getD(Position3d.class, e.getTargetVertex());
		double l = Rn.euclideanDistance(s, t);
		return getNewLength(e, l, u);
	}
	
	@Override
	public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
		return Edge.class.isAssignableFrom(nodeClass);
	}

	@Override
	public double getPriority() {
		return 10;
	}

}
