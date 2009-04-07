package de.varylab.discreteconformal.heds.util;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;

import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.TreeSet;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;

public class SpanningTreeUtility {

	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<E> getSpanningTree(Set<E> edges, E root) {
		Set<E> r = new TreeSet<E>(new NodeComparator<E>());
		if (root.getHalfEdgeDataStructure().numVertices() <= 1) {
			return r;
		}
		r.add(root);
		r.add(root.getOppositeEdge());
		Queue<E> growQ = new LinkedList<E>();
		growQ.offer(root);
		growQ.offer(root.getOppositeEdge());
		while (!growQ.isEmpty()) {
			E e = growQ.poll();
			Set<E> n = growTree(edges, r, e);
			for (E ne : n) {
				growQ.offer(ne);
			}
		}
		return r;
	}
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	>  Set<E> growTree(Set<E> edges, Set<E> tree, E branch) {
		Set<E> r = new TreeSet<E>(new NodeComparator<E>());
		V t = branch.getTargetVertex();
		List<E> star = incomingEdges(t);
		for (E e : star) {
			if (!edges.contains(e)) {
				continue;
			}
			Set<E> nextStar = new TreeSet<E>(new NodeComparator<E>());
			nextStar.addAll(incomingEdges(e.getStartVertex()));
			nextStar.retainAll(tree);
			if (nextStar.size() != 0) {
				continue;
			}
			tree.add(e);
			tree.add(e.getOppositeEdge());
			r.add(e);
			r.add(e.getOppositeEdge());
		}
		return r;
	}
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	>  Set<E> growTreeDual(Set<E> edges, Set<E> tree, E branch) {
		Set<E> r = new TreeSet<E>(new NodeComparator<E>());
		F t = branch.getLeftFace();
		List<E> star = HalfEdgeUtils.boundaryEdges(t);
		for (E e : star) {
			if (!edges.contains(e)) {
				continue;
			}
			F next = e.getRightFace();
			if (next == null) {
				continue;
			}
			Set<E> nextStar = new TreeSet<E>(new NodeComparator<E>());
			nextStar.addAll(HalfEdgeUtils.boundaryEdges(next));
			nextStar.retainAll(tree);
			if (nextStar.size() != 0) {
				continue;
			}
			tree.add(e);
			tree.add(e.getOppositeEdge());
			r.add(e);
			r.add(e.getOppositeEdge());
		}
		return r;
	}
	
	
	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<E> getDualSpanningTree(Set<E> edges, E root) {
		Set<E> r = new TreeSet<E>(new NodeComparator<E>());
		if (root.getHalfEdgeDataStructure().numFaces() <= 1) {
			return r;
		}
		r.add(root);
		r.add(root.getOppositeEdge());
		Queue<E> growQ = new LinkedList<E>();
		growQ.offer(root);
		growQ.offer(root.getOppositeEdge());
		while (!growQ.isEmpty()) {
			E e = growQ.poll();
			Set<E> n = growTreeDual(edges, r, e);
			for (E ne : n) {
				growQ.offer(ne);
			}
		}
		return r;
	}
	
}
