package de.varylab.discreteconformal.heds.util;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.jtem.halfedge.util.HalfEdgeUtils.isInteriorEdge;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.Stack;
import java.util.Vector;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;


public class Search {

	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Vector<E> bFS(V start, V end, boolean avoidBorder) throws NoSuchElementException{
		HashSet<V> endSet = new HashSet<V>();
		endSet.add(end);
		return bFS(start, endSet, avoidBorder);
	}

	/**
	 * Breadth first seach to the first hit in endPoints
	 * @param start
	 * @param end
	 * @param graph
	 * @return
	 * @throws NoSuchElementException
	 */
	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Vector<E> bFS(V start, Set<V> endPoints, boolean avoidBorder) throws NoSuchElementException{
		HashMap<V, Stack<E>> pathMap = new HashMap<V, Stack<E>>();
		LinkedList<V> queue = new LinkedList<V>();
		HashSet<V> visited = new HashSet<V>();
		V actVertex = start;
		queue.add(start);
		pathMap.put(start, new Stack<E>());
		while (!queue.isEmpty()){
			actVertex = queue.poll();
			Stack<E> path = pathMap.get(actVertex);
			List<E> star = incomingEdges(actVertex);
			for (E e : star){
				if (!isInteriorEdge(e) && avoidBorder)
					continue;
				E pathEdge = e.getOppositeEdge();
				V v = pathEdge.getTargetVertex();
				Stack<E> newPath = new Stack<E>();
				newPath.addAll(path);
				newPath.push(pathEdge);
				pathMap.put(v, newPath);
				if (!visited.contains(v)){
					visited.add(v);
					if (endPoints.contains(v))
						return newPath;
					else
						queue.offer(v);
				}
			}
		}
		throw new NoSuchElementException();
	}

	/**
	 * Depth first search - Untested!!
	 * @param start
	 * @param end
	 * @param graph
	 * @return
	 * @throws NoSuchElementException
	 */
	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Vector<E> dFS(V start, V end, HalfEdgeDataStructure<V, E, F> graph) throws NoSuchElementException{
		Stack<E> path = new Stack<E>();
		HashSet<V> visited = new HashSet<V>();
		visited.add(start);
		if (!advanceDFS(start, end, path, visited))
			throw new NoSuchElementException();
		else 
			return path;
	}

	private static 
	<
		V extends Vertex<V, E, ?>,
		E extends Edge<V, E, ?>
	> boolean advanceDFS(V start, V end, Stack<E> path, HashSet<V> visited) {
		if (start == end)
			return true;
		List<E> star = incomingEdges(start);
		for (E e : star){
			if (!visited.contains(e.getTargetVertex())){
				visited.add(e.getTargetVertex());
				path.push(e);
				if (advanceDFS(e.getTargetVertex(), end, path, visited))
					return true;
				path.pop();
			}
		}
		return false;
	}

}
