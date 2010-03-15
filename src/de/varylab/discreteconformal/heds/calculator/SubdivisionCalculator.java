package de.varylab.discreteconformal.heds.calculator;

import java.util.List;

import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.algorithm.calculator.EdgeAverageCalculator;
import de.jtem.halfedgetools.algorithm.calculator.FaceBarycenterCalculator;
import de.jtem.halfedgetools.algorithm.calculator.VertexPositionCalculator;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

public class SubdivisionCalculator implements VertexPositionCalculator, EdgeAverageCalculator, FaceBarycenterCalculator {

	private double
		alpha = 0.5;
	
	@Override
	public double getPriority() {
		return 0;
	}
	

	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double[] get(V v) {
		CoVertex cv = (CoVertex)v;
		return cv.getPosition().get();
	}

	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> void set(V v, double[] c) {
		CoVertex cv = (CoVertex)v;
		cv.getPosition().set(c);
	}
	
	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double[] get(E e) {
		CoEdge je = (CoEdge)e;
		double[] s = je.getStartVertex().getPosition().get();
		double[] t = je.getTargetVertex().getPosition().get();
		return Rn.linearCombination(null, alpha, t, 1 - alpha, s);
	}
	
	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double[] get(F f) {
		double[] pos = new double[3];
		List<E> b = HalfEdgeUtils.boundaryEdges(f);
		for (E e : b) {
			CoVertex jv = (CoVertex)e.getTargetVertex();
			Rn.add(pos, pos, jv.getPosition().get());
		}
		return Rn.times(pos, 1.0 / b.size(), pos);
	}

	@Override
	public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
		boolean result = false;
		result |= CoEdge.class.isAssignableFrom(nodeClass);
		result |= CoFace.class.isAssignableFrom(nodeClass);
		return result;
	}

	@Override
	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}

	@Override
	public void setIgnore(boolean ignore) {
	}

}
