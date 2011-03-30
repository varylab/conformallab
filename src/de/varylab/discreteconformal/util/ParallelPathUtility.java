package de.varylab.discreteconformal.util;

import java.util.List;
import java.util.Set;
import java.util.Vector;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.util.HalfEdgeUtilsExtra;

public class ParallelPathUtility {

	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> List<List<E>> getLeftParallelPaths(HalfEdgeDataStructure<V,E,F> hds, List<List<E>> paths) {
		List<List<E>> parallel= new Vector<List<E>>();
		for (List<E> path: paths) {
			parallel.add(getLeftParallelPath(hds, path));
		}
		return parallel;
	}
	
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> List<E> getLeftParallelPath(
			HalfEdgeDataStructure<V, E, F> hds, List<E> path) {
		
		Set<F> dualpathvertices = EdgeUtility.getDualVertexSet(DualityUtility
				.getDualPath(hds, path));
		List<E> boundary = getBoundary(dualpathvertices);
		List<E> inverseBoundary = invertPath(boundary);
		List<E> sum = addPaths(path, inverseBoundary);
		// TODO: delete the remaining triangle boundaries (spikes where deleted
		// from the dual path) 
		return sum;
	}
	
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> List<E> getBoundary(Set<F> faces){
		List<F> facelist= new Vector<F>();
		for (F f:faces) {
			facelist.add(f);
		}
		return getBoundary(facelist);
	}
	
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> List<E> getBoundary(List<F> faces){
		List<E> boundary= new Vector<E>();
		for (F f:faces) {
			List<E> b = HalfEdgeUtilsExtra.getBoundary(f);
			for (E e: b)
				boundary.add(e);
		}
		EdgeUtility.removeEdgePairs(boundary);
		return boundary;
	}
	
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> List<E> addPaths(List<E> path1,List<E> path2){
		List<E> sum= new Vector<E>();
		for (E e:path1)
			sum.add(e);
		for (E e:path2)
			sum.add(e);
		EdgeUtility.removeEdgePairs(sum);
		return sum;
	}
	
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> List<E> invertPath(List<E> path) {
		List<E> inv= new Vector<E>();
		for (E e: path)
			inv.add(e.getOppositeEdge());
		return inv;
	}
	
}
