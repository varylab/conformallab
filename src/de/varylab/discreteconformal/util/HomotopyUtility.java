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

public class HomotopyUtility {

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<Set<E>> getGeneratorPaths(V start, WeightAdapter<E> wa) {
		List<Set<E>> result = new ArrayList<Set<E>>();
		
		HalfEdgeDataStructure<V, E, F> hds = start.getHalfEdgeDataStructure();
		Set<E> edgeSet = new TreeSet<E>(new NodeIndexComparator<E>());
		edgeSet.addAll(hds.getEdges());
		Set<V> vSet = new TreeSet<V>(new NodeIndexComparator<V>());
		vSet.addAll(hds.getVertices());
		Set<E> tree1 = Search.getAllShortestPathsTree(start, vSet, wa, new HashSet<V>());
		edgeSet.removeAll(tree1);
		Set<E> tree2 = SpanningTreeUtility.getDualSpanningTree(edgeSet, edgeSet.iterator().next());
		edgeSet.removeAll(tree2);

		for (E bridge : edgeSet) {
			if (bridge.isPositive()) {
				Set<E> cycle = null;
				if (bridge.getStartVertex() == bridge.getTargetVertex()) {
					cycle = singleton(bridge);
				} else {
					cycle = HomologyUtility.findCycle(tree1, bridge);
				}
//				Set<V> vCycle = PathUtility.getVerticesOnPath(cycle);
//				if (!vCycle.contains(start)) { 
//					// cycle does not visit root, we have to add the connection manually
//					List<E> connectPath = Search.bFS(tree1, start, vCycle, true, null);
//					cycle.addAll(PathUtility.getFullPath(connectPath));
//				}
				result.add(cycle);
			}
		}
		return result;
	}
	
	
}
