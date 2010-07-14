package de.varylab.discreteconformal.util;

import static java.util.Collections.singleton;

import java.util.ArrayList;
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


	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<E> findCycle(Set<E> set, E bridge) {
		List<E> path = Search.bFS(set, bridge.getTargetVertex(), bridge.getStartVertex());
		Set<E> r = new TreeSet<E>(new NodeComparator<E>());
		r.addAll(path);
		r.add(bridge);
		return r;
	}
	
	
	public static
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<Set<E>> getGeneratorPaths(V root, WeightAdapter<E> wa) {
		List<Set<E>> result = new ArrayList<Set<E>>();
		
		HalfEdgeDataStructure<V, E, F> hds = root.getHalfEdgeDataStructure();
		Set<E> edgeSet = new TreeSet<E>(new NodeComparator<E>());
		edgeSet.addAll(hds.getEdges());
		Set<V> vSet = new TreeSet<V>(new NodeComparator<V>());
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
	

	
}
