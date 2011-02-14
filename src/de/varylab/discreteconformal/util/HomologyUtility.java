package de.varylab.discreteconformal.util;

import static java.util.Collections.singleton;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

public class HomologyUtility {

	protected static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<E> findCycle(Set<E> set, E bridge) {
		List<E> path = Search.bFS(set, bridge.getTargetVertex(), bridge.getStartVertex());
		Set<E> r = new TreeSet<E>(new NodeIndexComparator<E>());
		r.addAll(path);
		r.add(bridge);
		return r;
	}
	
	protected static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<E> findDualCycle(Set<E> set, E bridge) {
		List<E> path = Search.dualBFS(set, bridge.getLeftFace(), bridge.getRightFace());
		path.add(bridge);
		return path;
	}
	
	/**
	 * Calculates a basis for the homotopy of the given surface.
	 * The resulting cycles are visiting the root vertex, possibly
	 * over connecting paths.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param root The root of the system of cycles.
	 * @param wa An edge weight. Cycles minimize length with respect to this weights
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<Set<E>> getGeneratorPaths(V root, WeightAdapter<E> wa) {
		List<Set<E>> result = new ArrayList<Set<E>>();
		
		HalfEdgeDataStructure<V, E, F> hds = root.getHalfEdgeDataStructure();
		Set<E> edgeSet = new TreeSet<E>(new NodeIndexComparator<E>());
		edgeSet.addAll(hds.getEdges());
		Set<V> vSet = new TreeSet<V>(new NodeIndexComparator<V>());
		vSet.addAll(hds.getVertices());
		Set<E> tree1 = Search.getAllShortestPathsTree(root, vSet, wa, new HashSet<V>());
		edgeSet.removeAll(tree1);
		Set<E> tree2 = SpanningTreeUtility.getDualSpanningTree(edgeSet, edgeSet.iterator().next());
		edgeSet.removeAll(tree2);

		for (E bridge : edgeSet) {
			if (bridge.isPositive()) {
				if (bridge.getStartVertex() == bridge.getTargetVertex()) {
					result.add(singleton(bridge));
				} else {
					result.add(findCycle(tree1, bridge));
				}
			}
		}
		return result;
	}
	
	
	/**
	 * TODO Dual cycles might be disconnected 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param root
	 * @param wa
	 * @return
	 * @see getGeneratorPaths
	 */
	public static
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<List<E>> getDualGeneratorPaths(V root, WeightAdapter<E> wa) {
		Set<List<E>> result = new HashSet<List<E>>();
		
		HalfEdgeDataStructure<V, E, F> hds = root.getHalfEdgeDataStructure();
		Set<E> edgeSet = new TreeSet<E>(new NodeIndexComparator<E>());
		edgeSet.addAll(hds.getEdges());
		Set<V> vSet = new TreeSet<V>(new NodeIndexComparator<V>());
		vSet.addAll(hds.getVertices());
		Set<E> tree1 = Search.getAllShortestPathsTree(root, vSet, wa, new HashSet<V>());
		edgeSet.removeAll(tree1);
		Set<E> tree2 = SpanningTreeUtility.getDualSpanningTree(edgeSet, edgeSet.iterator().next());
		edgeSet.removeAll(tree2);

		for (E bridge : edgeSet) {
			if (bridge.isPositive()) {
				if (bridge.getLeftFace() == bridge.getRightFace()) {
					result.add(Collections.singletonList(bridge));
				} else {
					result.add(findDualCycle(tree2, bridge));
				}
			}
		}
		return result;
	}
}
