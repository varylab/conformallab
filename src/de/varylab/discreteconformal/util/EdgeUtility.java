package de.varylab.discreteconformal.util;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.EdgeIndex;


/**
 * Utility class for edges and cycles.
 * 
 * @author knoeppel
 * 
 */
public class EdgeUtility {

	/**
	 * There are several possibilities how an edge can be connected to the given
	 * cycle. If the surface is cut along the cycle, we can regard the cycle as
	 * two cycles actually - a left and a right one.
	 */
	public static enum EdgeStatus {
		endsAtLeftCycle,
		startsAtLeftCycle,
		endsAtRightCycle,
		startsAtRightCycle,
		liesOnCycle,
		liesOnInverseCycle,
		noConnection
	}

	/**
	 * Returns the intersection number of two oriented primal cycles.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param cycle1
	 * @param cycle2
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> int getIntersectionNumberOfPrimalCycles(List<E> cycle1, List<E> cycle2) {
		
		// get the number of positive intersections of the first cycle with the
		// second
		int numPositiveIntersections = countPrimalEdgesWithStatus(cycle1, cycle2,
				EdgeStatus.startsAtLeftCycle);
		// get the number of negative intersections of the first cycle with the
		// second
		int numNegativeIntersections = countPrimalEdgesWithStatus(cycle1, cycle2,
				EdgeStatus.endsAtLeftCycle);
		return numPositiveIntersections - numNegativeIntersections;
	}
	
	/**
	 * Returns the intersection number of two oriented dual cycles.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param cycle1
	 * @param cycle2
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> int getIntersectionNumberOfDualCycles(List<E> cycle1, List<E> cycle2) {
		
		// get the number of positive intersections of the first cycle with the
		// second
		int numPositiveIntersections = countDualEdgesWithStatus(cycle1, cycle2,
				EdgeStatus.startsAtLeftCycle);
		// get the number of negative intersections of the first cycle with the
		// second
		int numNegativeIntersections = countDualEdgesWithStatus(cycle1, cycle2,
				EdgeStatus.endsAtLeftCycle);
		return numPositiveIntersections - numNegativeIntersections;
	}

	/**
	 * Count the edges in a given list having a specified status with respect to
	 * a given primal cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param cycle
	 * @param list
	 * @param status
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> int countPrimalEdgesWithStatus(
				List<E> cycle, List<E> list, EdgeStatus status) {
		
		int count = 0;
		Set<V> vertexSet= getPrimalVertexSet(cycle);
		for(E e:list){
			if(getPrimalEdgeStatus(e, cycle,vertexSet)== status)
				count++;
		}
		return count;
	}
	
	/**
	 * Count the edges in a given list having a specified status with respect to
	 * a given dual cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param cycle
	 * @param list
	 * @param status
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> int countDualEdgesWithStatus(
				List<E> cycle, List<E> list, EdgeStatus status) {
		
		int count = 0;
		Set<F> vertexSet= getDualVertexSet(cycle);
		for(E e:list){
			if(getDualEdgeStatus(e, cycle,vertexSet)== status)
				count++;
		}
		return count;
	}

	/**
	 * Method gets a list of cycles and writes them into the columns of a matrix.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param adapters
	 * @param delaunay
	 * @param cycles
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix2D cyclesToMatrix(
			AdapterSet adapters, HalfEdgeDataStructure<V, E, F> delaunay,
			List<List<E>> cycles) {
		
		// initialize two matrices to encode the cycles
		DoubleMatrix2D A = new SparseDoubleMatrix2D(delaunay.numEdges() / 2,
				cycles.size());

		// needed to fill the matrices
		Integer currEdgeId;
		double currval;

		// fill the matrix of a-cycles
		for (int i = 0; i < cycles.size(); i++) {
			for (E e : cycles.get(i)) {
				currEdgeId = adapters.get(EdgeIndex.class, e, Integer.class);
				currval = A.get(currEdgeId, i);
				if (e.isPositive())
					currval++;
				else
					currval--;
				A.set(currEdgeId, i, currval);
			}
		}
		// print(A, 0);

		return A;
	}
	
	/**
	 * Method gets a matrix of cycles and returns them as a list.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param adapters
	 * @param delaunay
	 * @param cycles
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<List<E>> matrixToCycles(
			AdapterSet adapters, HalfEdgeDataStructure<V, E, F> delaunay,
			DoubleMatrix2D cycles) {
		
		// initialize two matrices to encode the cycles
		List<List<E>> A = new Vector<List<E>>();
		
		int n= cycles.columns();
		for (int i = 0; i < n; i++) {
			A.add(new Vector<E>());
		}

		// needed to fill the matrices
		Integer posEdgeId;
		double currval;

		// iterate over positive edges to build up the new cycles
		for (E e : delaunay.getPositiveEdges()) {
			// get positive edge index
			posEdgeId = adapters.get(EdgeIndex.class, e, Integer.class);
			// add edge with this index to the cycles
			for (int i = 0; i < n; i++) {
				currval = cycles.get(posEdgeId, i);
				// add the positive edge if the weight is positive else the
				// opposite edge
				E currEdge = (currval > 0) ? e : e.getOppositeEdge();
				// as many times as the weight says
				for (int count = 0; count < Math.abs(currval); count++)
					A.get(i).add(currEdge);
			}
		}

		return A;
	}
	
	/**
	 * Rotate e around v clockwise
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param e
	 * @param v
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> E getNextPrimalEdgeClockwise(E e, V v) {
		
		if (e.getStartVertex() == v) {
			return e.getOppositeEdge().getNextEdge();
		}
		if (e.getTargetVertex() == v) {
			return e.getNextEdge().getOppositeEdge();
		}
		throw new IllegalArgumentException(
				"Edge does not contain vertex in getNextEdgeClockwise()");
	}
		
	/**
	 * Rotate e around f clockwise
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param e
	 * @param f
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> E getNextDualEdgeClockwise(E e, F f) {
		
		if (e.getLeftFace() == f) {
			return e.getPreviousEdge();
		}
		if (e.getRightFace() == f) {
			return e.getOppositeEdge().getPreviousEdge().getOppositeEdge();
		}
		throw new IllegalArgumentException(
				"Edge does not contain vertex in getNextEdgeClockwise()");
	}

	/**
	 * Returns the edge's status, see EdgeStatus description.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param e
	 * @param edgeCycle
	 * @param vertexCycle
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> EdgeStatus getPrimalEdgeStatus(E e, List<E> edgeCycle, Set<V> vertexCycle) {
		
		if (edgeCycle.contains(e))
			return EdgeStatus.liesOnCycle;
		if (edgeCycle.contains(e.getOppositeEdge()))
			return EdgeStatus.liesOnInverseCycle;
		boolean outpointing = vertexCycle.contains(e.getStartVertex());
		V v = outpointing ? e.getStartVertex() : e.getTargetVertex();
		E curr = getNextPrimalEdgeClockwise(e, v);
		while (curr != e) {
			if (edgeCycle.contains(curr)) {
				if (outpointing)
					return EdgeStatus.startsAtLeftCycle;
				else
					return EdgeStatus.endsAtRightCycle;
			}
			if (edgeCycle.contains(curr.getOppositeEdge())) {
				if (outpointing)
					return EdgeStatus.startsAtRightCycle;
				else
					return EdgeStatus.endsAtLeftCycle;
			}
			curr = getNextPrimalEdgeClockwise(curr, v);
		}

		return EdgeStatus.noConnection;
	}
	
	/**
	 * Returns the dual edge's status, see EdgeStatus description.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param e
	 * @param edgeCycle
	 * @param dualVertexCycle
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> EdgeStatus getDualEdgeStatus(E e, List<E> edgeCycle, Set<F> dualVertexCycle) {
		
		if (edgeCycle.contains(e))
			return EdgeStatus.liesOnCycle;
		if (edgeCycle.contains(e.getOppositeEdge()))
			return EdgeStatus.liesOnInverseCycle;
		boolean outpointing = dualVertexCycle.contains(e.getRightFace());
		F f = outpointing ? e.getRightFace() : e.getLeftFace();
		E curr = getNextDualEdgeClockwise(e, f);
		while (curr != e) {
			if (edgeCycle.contains(curr)) {
				if (outpointing)
					return EdgeStatus.startsAtRightCycle;
				else
					return EdgeStatus.endsAtLeftCycle;
			}
			if (edgeCycle.contains(curr.getOppositeEdge())) {
				if (outpointing)
					return EdgeStatus.startsAtLeftCycle;
				else
					return EdgeStatus.endsAtRightCycle;
			}
			curr = getNextDualEdgeClockwise(curr, f);
		}

		return EdgeStatus.noConnection;
	}
	
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> void removeEdgePairs(
				List<E> edgeCycle) {
		List<E> pairs= new Vector<E>();
		for(E e: edgeCycle){
			if(edgeCycle.contains(e.getOppositeEdge())){
				pairs.add(e); pairs.add(e.getOppositeEdge());
			}
		}
		for (E e: pairs) {
			edgeCycle.remove(e);
		}
	}
	
	/**
	 * Returns the indices of the vertices contained in the cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param cycle
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<V> getPrimalVertexSet(List<E> cycle) {
		
		Set<V> vertexSet = new HashSet<V>(cycle.size());
		for (E e:cycle) {
			vertexSet.add(e.getStartVertex());
		}
		return vertexSet;
	}
	
	/**
	 * Returns the indices of the faces contained in the dual cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param cycle
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<F> getDualVertexSet(List<E> cycle) {
		
		Set<F> vertexSet = new HashSet<F>(cycle.size());
		for (E e:cycle) {
			vertexSet.add(e.getRightFace());
		}
		return vertexSet;
	}

}
