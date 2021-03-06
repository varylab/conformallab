package de.varylab.discreteconformal.util;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryEdge;
import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryVertex;
import static de.jtem.halfedge.util.HalfEdgeUtils.isInteriorEdge;
import static de.jtem.halfedge.util.HalfEdgeUtils.isManifoldVertex;
import static java.lang.Math.abs;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition2d;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

public class CuttingUtility {

	private static Logger
		log = Logger.getLogger(CuttingUtility.class.getName());
	
	public static class CuttingInfo <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> {
		public V
			cutRoot = null;
		public Map<E, E>
			edgeCutMap = new HashMap<E, E>();
		public Map<V, V>
			vertexCopyMap = new HashMap<V, V>();
		
		private Set<Set<E>>
			_paths = new HashSet<Set<E>>();
		private Map<Set<E>, Set<E>>
			_pathCutMap = new HashMap<Set<E>, Set<E>>();
		
		public Set<V> getCopies(V v) {
			Set<V> copies = new TreeSet<V>(new NodeIndexComparator<V>());
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
			Set<V> branches = new TreeSet<V>(new NodeIndexComparator<V>());
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
	> void glueAll(HDS hds, CuttingInfo<V, E, F> cut) {
		for (E e1 : cut.edgeCutMap.keySet()) {
			E e2 = cut.edgeCutMap.get(e1);
			if (e1.isValid() && e2.isValid() && isBoundaryEdge(e1) && isBoundaryEdge(e2)) {
				glue(e1, e2);
			}
		}
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void glue(E e1, E e2) {
		if (!isBoundaryEdge(e1) || !isBoundaryEdge(e2)) {
			throw new RuntimeException("Can only glue boundary edges");
		}
		if (e1.getLeftFace() == null) {
			e1 = e1.getOppositeEdge();
		}
		if (e2.getLeftFace() == null) {
			e2 = e2.getOppositeEdge();
		}
		if (e1.getStartVertex() != e2.getTargetVertex()) {
			StitchingUtility.stitch(e1.getStartVertex(), e2.getTargetVertex());
		}
		if (e1.getTargetVertex() != e2.getStartVertex()) {
			StitchingUtility.stitch(e1.getTargetVertex(), e2.getStartVertex());
		}
	}
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> CuttingInfo<V, E, F> cutTorusToDisk(HDS hds, V root, WeightAdapter<E> wa) {
		if (HalfEdgeUtils.getGenus(hds) != 1) throw new RuntimeException("Invalid genus in cutTorusToDisk()");
		CuttingInfo<V, E, F> context = new CuttingInfo<V, E, F>();
		List<Set<E>> paths = HomotopyUtility.getGeneratorPaths(root, wa);
		Set<E> path0 = paths.get(0);
		for (Set<E> path : paths) {
			if (path.size() < path0.size()) {
				path0 = path;
			}
		}
		context._paths.add(path0);
		cutAlongPath(path0, context);

		List<E> path1 = null;
		V cutRoot = null;
		Collection<V> boundary = HalfEdgeUtils.boundaryVertices(hds);
		List<V> checkList = new LinkedList<V>();
		checkList.add(root); // we want the root to be the cut root if possible
		checkList.addAll(context.vertexCopyMap.keySet());
		for (V bv : checkList) {
			V cbv = context.vertexCopyMap.get(bv);
			Set<V> avoidSet = new HashSet<V>(boundary);
			avoidSet.remove(bv);
			avoidSet.remove(cbv);
			try {
				path1 = Search.bFS(bv, cbv, avoidSet);
			} catch (Exception e) {
				continue;
			}
			cutRoot = bv;
			break;
		}
		assert path1 != null;
		assert cutRoot != null;
		
		context.cutRoot = cutRoot;
		context._paths.add(new HashSet<E>(path1));
		cutAlongPath(path1, context);
		for (Set<E> path : context._paths) {
			Set<E> coPath = new TreeSet<E>(new NodeIndexComparator<E>());
			context._pathCutMap.put(path, coPath);
			for (E e : path) {
				E coE = context.edgeCutMap.get(e);
				if (coE != null) {
					coPath.add(coE);
				}
			}
		}
		return context;
	}
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	>  void cutAlongPath(Collection<E> path, CuttingInfo<V, E, F> context) {
		
		for (E e : path) {
			if(context.cutRoot == null) {
				context.cutRoot = e.getStartVertex();
			}
			if (isInteriorEdge(e)) {
				E eop = e.getOppositeEdge();
				context.edgeCutMap.put(e, eop);
				context.edgeCutMap.put(eop, e);
				Map<V, V> vMap = cutAtEdge(e);
				context.edgeCutMap.put(e.getOppositeEdge(), eop.getOppositeEdge());
				context.edgeCutMap.put(eop.getOppositeEdge(), e.getOppositeEdge());
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
	}
	
	
	
	/**
	 * Cuts a manifold of genus g along 2*g generator paths.
	 * In general the paths share edges.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param <HDS>
	 * @param hds
	 * @param root
	 * @param wa
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> CuttingInfo<V, E, F> cutManifoldToDisk(HDS hds, V root, WeightAdapter<E> wa) {
		int genus = HalfEdgeUtils.getGenus(hds);
		CuttingInfo<V, E, F> context = new CuttingInfo<V, E, F>();
		context.cutRoot = root;
		if (genus < 1) {
			return context;
		}
		log.info("cutting manifold to disk...");
		log.info("cut root is " + root);
		for(Set<E> path : HomotopyUtility.getGeneratorPaths(root, wa)) {
			log.info("path: " + path);
			context._paths.add(path);
		}
		Set<E> masterPath = new TreeSet<E>(new NodeIndexComparator<E>());
		for (Set<E> path : context._paths) {
			masterPath.addAll(path);
		}
		log.info("cutting along edges: " + masterPath);
		cutAtEdges(context, masterPath);
		for (Set<E> path : context._paths) {
			Set<E> coPath = new TreeSet<E>(new NodeIndexComparator<E>());
			context._pathCutMap.put(path, coPath);
			for (E e : path) {
				E coE = context.edgeCutMap.get(e);
				if (coE != null) {
					coPath.add(coE);
				}
			}
		}
		return context;
	}


	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> void cutAtEdges(CuttingInfo<V, E, F> context, Set<E> cutEdges) {
		for (E e : cutEdges) { 
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
		new1.copyData(edge);
		new2.copyData(opp);
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
			if (b == null) {
				throw new RuntimeException("Could not find an edge at the boundary in cutAtEdge().");
			}
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
			if (b == null) {
				throw new RuntimeException("Could not find an edge at the boundary in cutAtEdge().");
			}
			E b2 = b.getNextEdge();
			List<E> newTargetEdges = new LinkedList<E>();
			E actEdge = new1;
			do {
				newTargetEdges.add(actEdge);
				actEdge = actEdge.getOppositeEdge().getPreviousEdge();
			} while (actEdge != b2.getOppositeEdge());
			newTargetEdges.add(b2.getOppositeEdge());
			
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
			if (be1 == null) {
				throw new RuntimeException("Could not find an edge at the boundary in cutLoopEdge().");
			}
			
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

	public static void cutToSimplyConnected(CoHDS hds, CoVertex cutRoot, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		List<List<CoEdge>> bc = HalfEdgeUtils.boundaryComponents(hds);
		if (cutRoot != null && bc.size() > 1) {
			for (List<CoEdge> path : bc) {
				Set<CoVertex> cycle = PathUtility.getVerticesOnPath(path);
				List<CoEdge> cutPath = Search.getShortestPath(cutRoot, cycle, null);
				cutAlongPath(cutPath, cutInfo);
			}
			cutInfo.cutRoot = cutRoot;
		} else {
			cutInfo.cutRoot = bc.get(0).get(0).getStartVertex();
			while (bc.size() > 1) {
				Set<CoVertex> cycle1 = PathUtility.getVerticesOnPath(bc.get(0));
				List<CoEdge> allOtherCycles = new LinkedList<CoEdge>();
				for (int i = 1; i < bc.size(); i++) {
					allOtherCycles.addAll(bc.get(i));
				}
				Set<CoVertex> otherVertices = PathUtility.getVerticesOnPath(allOtherCycles);
				
				Set<CoVertex> inter = new TreeSet<CoVertex>(new NodeIndexComparator<CoVertex>());
				inter.addAll(cycle1);
				inter.retainAll(otherVertices);
				if (!inter.isEmpty()) {
					throw new RuntimeException("paths cannot intersect!");
				}
				List<CoEdge> cutPath = new LinkedList<CoEdge>();
				for (CoVertex s : cycle1) {
					List<CoEdge> checkPath = Search.getShortestPath(s, otherVertices, null);
					if (checkPath.size() < cutPath.size() || cutPath.isEmpty()) {
						cutPath = checkPath;
					}
				}
				if (cutPath.isEmpty()) {
					throw new RuntimeException("no path between cycles found!");
				}
				cutAlongPath(cutPath, cutInfo);
				bc = HalfEdgeUtils.boundaryComponents(hds);
			}
		}
	}

	
	/**
	 * Checks if the given p is inside the texture face of f.
	 * The texture points of f are assumed to build a convex polygon.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param <HDS>
	 * @param p
	 * @param f
	 * @param a
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> boolean isInConvexTextureFace(double[] p, F f, AdapterSet a) {
		List<E> bList = HalfEdgeUtils.boundaryEdges(f);
		boolean sideFlag = false;
		boolean firstEdge = true;
		for (E e : bList) {
			double[] sp = a.getD(TexturePosition2d.class, e.getStartVertex());
			double[] tp = a.getD(TexturePosition2d.class, e.getTargetVertex());
			double[] v1 = Rn.subtract(null, tp, sp);
			double[] v2 = Rn.subtract(null, p, sp);
			double[] v2r = {-v2[1], v2[0], 0};
			Rn.normalize(v1, v1);
			Rn.normalize(v2r, v2r);
			double dot = Rn.innerProduct(v1, v2r);
			if(abs(dot) <= 1E-4) {
				continue;
			}
			if (firstEdge) {
				sideFlag = dot > 0;
				firstEdge = false;
			} else {
				if (dot > 0 != sideFlag) {
					return false;
				}
			}
		}
		return true;
	}	
	
}
