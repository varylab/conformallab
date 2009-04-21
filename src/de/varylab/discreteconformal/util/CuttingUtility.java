package de.varylab.discreteconformal.util;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryVertex;
import static de.jtem.halfedge.util.HalfEdgeUtils.isInteriorEdge;
import static de.jtem.halfedge.util.HalfEdgeUtils.isManifoldVertex;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

public class CuttingUtility {

	
	public static class CuttingInfo <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> {
		
		public V
			cutRoot = null;
		public Map<E, E>
			edgeCutMap = new HashMap<E, E>();
		public List<Set<E>>
			paths = new ArrayList<Set<E>>();
		public Map<Set<E>, Set<E>>
			pathCutMap = new HashMap<Set<E>, Set<E>>();
		public Map<V, V>
			vertexCopyMap = new HashMap<V, V>();
		
		public Set<V> getCopies(V v) {
			Set<V> copies = new TreeSet<V>(new NodeComparator<V>());
			for (V cV : vertexCopyMap.keySet()) {
				V tmpV = vertexCopyMap.get(cV);
				while (tmpV != null) {
					if (tmpV == v) {
						copies.add(cV);
					}
					tmpV = vertexCopyMap.get(tmpV);
				}
			}
			V tmpV = vertexCopyMap.get(v);
			while (tmpV != null) {
				copies.add(tmpV);
				tmpV = vertexCopyMap.get(tmpV);
			}		
			copies.add(v);
			return copies;
		}
		
		public Set<V> getBranchSet() {
			Set<V> branches = new TreeSet<V>(new NodeComparator<V>());
			for (V v : vertexCopyMap.keySet()) {
				V copy = vertexCopyMap.get(v);
				V branchCopy = vertexCopyMap.get(copy);
				if (branchCopy != null) {
					branches.add(v);
					branches.add(copy);
					branches.add(branchCopy);
				}
			}
			return branches;
		}
	
	}
	
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> CuttingInfo<V, E, F> cutManifoldToDisk(HDS hds, V root, WeightAdapter<E> wa) {
		CuttingInfo<V, E, F> context = new CuttingInfo<V, E, F>();
		context.cutRoot = root;
		context.paths = HomologyUtility.getGeneratorPaths(root, wa);
		Set<E> masterPath = new TreeSet<E>(new NodeComparator<E>());
		for (Set<E> path : context.paths) {
			masterPath.addAll(path);
		}
		for (E e : masterPath) { 
			if (isInteriorEdge(e)) {
				context.edgeCutMap.put(e, e.getOppositeEdge());
				context.edgeCutMap.put(e.getOppositeEdge(), e);
				Map<V, V> vMap = cutAtEdge(e);
				for (V v : vMap.keySet()) {
					V copy = vMap.get(v);
					if (context.vertexCopyMap.keySet().contains(v)) {
						V oldCopy = context.vertexCopyMap.get(v);
						context.vertexCopyMap.put(copy, oldCopy);
					}
					context.vertexCopyMap.put(v, copy);
				}
			}
		}
		for (Set<E> path : context.paths) {
			Set<E> coPath = new TreeSet<E>(new NodeComparator<E>());
			context.pathCutMap.put(path, coPath);
			for (E e : path) {
				E coE = context.edgeCutMap.get(e);
				if (coE != null) {
					coPath.add(coE);
				}
			}
		}
		return context;
	}
	
	
	
	/**
	 * Cuts an edge and returns the mapping 
	 * oldVertex -> newVertex for split vertices
	 * @param edge
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Map<V, V> cutAtEdge(E edge) throws IllegalArgumentException{
		Map<V, V> result = new HashMap<V, V>();
		HalfEdgeDataStructure<V, E, F> hds = edge.getHalfEdgeDataStructure();
		V v1 = edge.getStartVertex();
		V v2 = edge.getTargetVertex();
		E opp = edge.getOppositeEdge();
		if (v1 == v2) {
			return cutLoopEdge(edge);
		}
		
		boolean splitV1 = isBoundaryVertex(v1);
		boolean splitV2 = isBoundaryVertex(v2);
		
		List<E> v1Star = incomingEdges(v1);
		List<E> v2Star = incomingEdges(v2);
		
		E new1 = hds.addNewEdge();
		E new2 = hds.addNewEdge();
		
		new1.linkOppositeEdge(opp);
		new2.linkOppositeEdge(edge);
		new1.setTargetVertex(v2);
		new2.setTargetVertex(v1);
		
		new1.linkNextEdge(new2);
		new2.linkNextEdge(new1);
		
		if (splitV1){
			E b = null;
			for (E e : v1Star) {
				if (e.getLeftFace() == null){
					b = e;
					break;
				}
			}
			assert (b != null);
			E b2 = b.getNextEdge();
			List<E> newTargetEdges = new LinkedList<E>();
			E actEdge = opp;
			do {
				newTargetEdges.add(actEdge);
				actEdge = actEdge.getNextEdge().getOppositeEdge();
			} while (actEdge != b);
			newTargetEdges.add(b);
			
			V newV = hds.addNewVertex();
			newV.copyData(v1);
			result.put(v1, newV);
			
			b.linkNextEdge(new1);
			new2.linkNextEdge(b2);
			
			for (E e : newTargetEdges) {
				e.setTargetVertex(newV);
			}
		}
		if (splitV2){
			E b = null;
			for (E e : v2Star) {
				if (e.getLeftFace() == null){
					b = e;
					break;
				}
			}
			assert (b != null);
			E b2 = b.getNextEdge();
			List<E> newTargetEdges = new LinkedList<E>();
			E actEdge = edge;
			do {
				newTargetEdges.add(actEdge);
				actEdge = actEdge.getNextEdge().getOppositeEdge();
			} while (actEdge != b);
			newTargetEdges.add(b);
			
			V newV = hds.addNewVertex();
			newV.copyData(v2);
			result.put(v2, newV);
			
			b.linkNextEdge(new2);
			new1.linkNextEdge(b2);
			
			for (E e : newTargetEdges) {
				e.setTargetVertex(newV);
			}
		}
		return result;
		
	}
	
	
	
	/**
	 * Cuts an edge which is a loop and returns the mapping 
	 * oldVertex -> newVertex for split vertices
	 * @param edge
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static 
	<
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Map<V, V> cutLoopEdge(E edge) throws IllegalArgumentException{
		Map<V, V> result = new HashMap<V, V>();
		HalfEdgeDataStructure<V, E, F> hds = edge.getHalfEdgeDataStructure();
		V oldVertex = edge.getStartVertex();
		V v2 = edge.getTargetVertex();
		if (oldVertex != v2) {
			throw new IllegalArgumentException("Start vertex != target vertex in cutLootEdge()");
		}
		
		boolean repairNonManifold = isBoundaryVertex(oldVertex);
		
		E opp = edge.getOppositeEdge();
		V newVertex = hds.addNewVertex();
		newVertex.copyData(oldVertex);
		result.put(oldVertex, newVertex);
		
		// get edges of the new vertex 
		List<E> newTargetEdges = new LinkedList<E>();
		E actEdge = edge;
		do {
			newTargetEdges.add(actEdge);
			actEdge = actEdge.getNextEdge().getOppositeEdge();
			// if this is no loop we have a dead lock here
		} while (actEdge != edge.getOppositeEdge()); 
		
		// create new edges
		E new1 = hds.addNewEdge();
		E new2 = hds.addNewEdge();
		
		new1.linkOppositeEdge(opp);
		new2.linkOppositeEdge(edge);
		new1.setTargetVertex(oldVertex);
		new2.setTargetVertex(newVertex);
		
		for (E e : newTargetEdges) {
			e.setTargetVertex(newVertex);
		}
		
		// link boundary loops
		new1.linkNextEdge(new1);
		new2.linkNextEdge(new2);

		// repair the created non manifold vertex
		if (repairNonManifold) {
			V nmv = null;
			if (isManifoldVertex(oldVertex)) {
				nmv = newVertex;
			} else {
				nmv = oldVertex;
			}
			assert isManifoldVertex(nmv == oldVertex ? newVertex : oldVertex);
			
			V v = hds.addNewVertex();
			v.copyData(oldVertex);
			result.put(oldVertex, v);
			
			E be1 = null;
			List<E> star = incomingEdges(nmv);
			for (E e : star) {
				if (e.getLeftFace() == null) {
					be1 = e;
					break;
				}
			}
			assert be1 != null;
			
			E be2 = be1.getNextEdge();
			actEdge = be1;
			do {
				actEdge.setTargetVertex(v);
				actEdge = actEdge.getOppositeEdge().getPreviousEdge();
			} while (actEdge.getLeftFace() != null);
			E be4 = actEdge;
			E be3 = be4.getNextEdge();
			
			be1.linkNextEdge(be3);
			be4.linkNextEdge(be2);
		}
		return result;
	}

	
}
