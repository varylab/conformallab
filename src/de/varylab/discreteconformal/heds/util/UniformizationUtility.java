package de.varylab.discreteconformal.heds.util;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.jtem.halfedge.util.HalfEdgeUtils.neighboringVertices;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.sqrt;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.functional.conformal.ConformalAdapters.Alpha;
import de.jtem.halfedge.functional.conformal.ConformalAdapters.Lambda;

public class UniformizationUtility {

	
	public static interface EdgeFlipAdapter<E extends Edge<?, E, ?>> {
		public void flip(E edge);
	}
	
	public static interface UAdapter<V extends Vertex<V, ?, ?>> {
		public double getU(V v);
	}
	

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void reduceToFundamentalPolygon(
		HDS hds, 
		V root,
		UAdapter<V> u,
		EdgeFlipAdapter<E> flip,
		Lambda<E> lambda,
		Alpha<E> alpha
	) {
		while (hds.numVertices() > 1) {
			boolean removed = false;
			for (V v : hds.getVertices()) {
				if (v == root) {
					continue;
				}
				int starSize = neighboringVertices(v).size();
				if (starSize == 3) {
					removeTrivalentVertex(v);
					removed = true;
				}
			}
			V v = root.getIncomingEdge().getStartVertex();
			if (!removed && hds.numVertices() != 1) {
				while (incomingEdges(v).size() != 3) {
					flip.flip(v.getIncomingEdge());
				}
			}
			removeTrivalentVertex(v);
		}
	}
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		U extends UAdapter<V>
	> double getNewLambda(
		E e,
		U u,
		Lambda<E> lambda,
		Alpha<E> alpha	
	) {
//		double a = e.getNextEdge().
//		double newLength = 
//		
//		
		return 0.0;
	}
	
	
	/**
	 * Calculate the edge length for the flat metric
	 * @param e
	 * @param u
	 * @return the new edge length
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		U extends UAdapter<V>
	> double getNewLength(
		E e, 
		U u,
		Lambda<E> lambda,
		Alpha<E> alpha
	) {
		V v1 = e.getStartVertex();
		V v2 = e.getTargetVertex();
		Double u1 = u.getU(v1); 
		Double u2 = u.getU(v2);
		Double la = lambda.getLambda(e);
		Double lambdaNew = la + u1 + u2;
		return 2 * arsinh( exp(lambdaNew / 2) );
	}
	
	
	private static double arsinh(double x) {
		double r = x + sqrt(x*x + 1);
		return log(r);
	}
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void removeTrivalentVertex(V v) {
		HalfEdgeDataStructure<V, E, F> hds = v.getHalfEdgeDataStructure();
		E in1 = v.getIncomingEdge();
		E out1 = in1.getOppositeEdge();
		E outer1 = in1.getPreviousEdge();
		E in2 = out1.getPreviousEdge();
		E out2 = in2.getOppositeEdge();
		E outer2 = in2.getPreviousEdge();
		E in3 = out2.getPreviousEdge();
		E out3 = in3.getOppositeEdge();
		E outer3 = in3.getPreviousEdge();
		
		F f1 = in1.getLeftFace();
		F f2 = in2.getLeftFace();
		F f3 = in3.getLeftFace();

		hds.removeVertex(v);
		hds.removeEdge(in1);
		hds.removeEdge(in2);
		hds.removeEdge(in3);
		hds.removeEdge(out1);
		hds.removeEdge(out2);
		hds.removeEdge(out3);
		hds.removeFace(f2);
		hds.removeFace(f3);
		
		outer1.setLeftFace(f1);
		outer2.setLeftFace(f1);
		outer3.setLeftFace(f1);
		
		outer1.linkNextEdge(outer2);
		outer2.linkNextEdge(outer3);
		outer3.linkNextEdge(outer1);
	}
	
	
	
}
