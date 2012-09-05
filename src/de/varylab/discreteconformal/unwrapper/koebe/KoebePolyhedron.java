package de.varylab.discreteconformal.unwrapper.koebe;

import static java.lang.Math.PI;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.generic.Position2d;
import de.jtem.halfedgetools.algorithm.subdivision.MedialGraphLinear;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.varylab.discreteconformal.unwrapper.circlepattern.CirclePatternLayout;
import de.varylab.discreteconformal.unwrapper.circlepattern.CirclePatternUtility;

/**
 * Calculates the koebe polyhedron
 * <p>
 * Copyright 2005 <a href="http://www.sechel.de">Stefan Sechelmann</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Stefan Sechelmann
 */
public class KoebePolyhedron {

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> KoebePolyhedronContext<V, E, F> contructKoebePolyhedron(HDS hds, AdapterSet a) {
		KoebePolyhedronContext<V, E, F> context = new KoebePolyhedronContext<V, E, F>();
		if (!HalfEdgeUtils.isValidSurface(hds)) {
			throw new RuntimeException("No valid surface in constructKoebePolyhedron()");
		}
		// medial graph
		HashMap<V, F> vertexFaceMap = new HashMap<V, F>();
		HashMap<F, F> faceFaceMap = new HashMap<F, F>();
		HashMap<E, E> edgeEdgeMap = new HashMap<E, E>();
		HashMap<E, V> edgeVertexMap = new HashMap<E, V>();
		HDS medial = instantiate(hds);
		MedialGraphLinear medialSubd = new MedialGraphLinear();
		medialSubd.execute(hds, medial, vertexFaceMap, edgeVertexMap, faceFaceMap, edgeEdgeMap, a);
		HalfEdgeUtils.fillAllHoles(medial);
		
		if (!HalfEdgeUtils.isValidSurface(medial)) {
			throw new RuntimeException("No valid medial graph could be contructed in constructKoebePolyhedron()");
		}
		
		// store old medial combinatorics
		HashMap<V, V> medialVertexMap = new HashMap<V, V>();
		HashMap<E, E> medialEdgeMap = new HashMap<E, E>();
		HashMap<F, F> medialFaceMap = new HashMap<F, F>();
		
		HalfEdgeDataStructure<V, E, F> medialCirclePattern = createCopy(medial, medialVertexMap, medialEdgeMap, medialFaceMap);
		
		// remove north pole to make a disk
		TopologyAlgorithms.removeVertex(medialCirclePattern.getVertex(0));
		
		// we want an orthogonal circle pattern
		Map<E, Double> thetaMap = new HashMap<E, Double>();
		for (E e : medialCirclePattern.getEdges()) {
			thetaMap.put(e, PI / 2);
		}
		Map<F, Double> phiMap = new HashMap<F, Double>();
		for (F f : medialCirclePattern.getFaces()) {
			phiMap.put(f, 2*PI);
		}
		
		if (!HalfEdgeUtils.isValidSurface(medialCirclePattern)) {
			throw new RuntimeException("No surface after altering medial graph in constructKoebePolyhedron()");
		}
		
		// optimization
		Map<F, Double> rhoMap = CirclePatternUtility.calculateCirclePatternRhos(medialCirclePattern, thetaMap, phiMap);
		CirclePatternLayout.execute(medialCirclePattern, rhoMap, thetaMap, a, Position2d.class, Position.class);
		
		// check rhos
		for (F f : medialCirclePattern.getFaces()){
			if (rhoMap.get(f) < -10 || rhoMap.get(f) > 10) {
				System.err.println("WARNING: some radii seem too big");
			}
		}
		
		normalizeBeforeProjection(medialCirclePattern, rhoMap, a);
		
		// projection faces
		HashSet<F> readyFaces = new HashSet<F>();
		for (E e : medial.getEdges()){
			E circlePatternEdge = medialEdgeMap.get(e);
			if (!circlePatternEdge.isValid()) {
				continue;
			}
			F f = e.getLeftFace();
			if (readyFaces.contains(f)) { 
				continue;
			}
			double[] peak = calculateConePeek(circlePatternEdge, rhoMap, a);
			a.set(Position.class, f, peak);
			readyFaces.add(f);
		}
		// projection vertices
		inverseStereographicProjection(medialCirclePattern, 1.0, a);
		
		// fill medial information
		V medialNorth = medial.getVertex(0);
		for (V v : medial.getVertices()){
			if (v == medialNorth) {
				continue;
			}
			V vertex = medialVertexMap.get(v);
			double[] p = a.getD(Position.class, vertex);
			a.set(Position.class, v, p);
		}
		a.set(Position.class, medialNorth, new double[] {0, 1, 0, 1});
		
		// fill into graph
		for (F face : hds.getFaces()){
			F f = faceFaceMap.get(face);
			double[] p = a.getD(Position.class, f);
			Pn.dehomogenize(p, p);
			a.set(Position.class, face, p);
		}
		for (V vertex : hds.getVertices()){
			F f = vertexFaceMap.get(vertex);
			double[] p = a.getD(Position.class, f);
			Pn.dehomogenize(p, p);
			a.set(Position.class, vertex, p);
		}
		
		// context
		context.edgeEdgeMap = edgeEdgeMap;
		context.edgeVertexMap = edgeVertexMap;
		context.faceFaceMap = faceFaceMap;
		context.vertexFaceMap = vertexFaceMap;
		context.circlePattern = medialCirclePattern;
		context.medial = medial;
		context.northPole = medialNorth;
		context.polyhedron = hds;
		return context;
	}


	@SuppressWarnings("unchecked")
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	>  HDS instantiate(HDS hds) {
		HDS newHDS = null;
		try {
			newHDS = (HDS)hds.getClass().newInstance();
		} catch (Exception e) {
			throw new RuntimeException("could not create an empty data structure from " + hds);
		}
		return newHDS;
	}
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		T extends HalfEdgeDataStructure<V, E, F>
	> HalfEdgeDataStructure<V, E, F> createCopy(T heds, HashMap<V, V> vertexMap, HashMap<E, E> edgeMap, HashMap<F, F> faceMap){
		// cloning nodes
		vertexMap.clear();
		edgeMap.clear();
		faceMap.clear();
		T copy = instantiate(heds);
		heds.createCombinatoriallyEquivalentCopy(copy);
		for (V v : copy.getVertices()) {
			vertexMap.put(heds.getVertex(v.getIndex()), v);
		}
		for (E e : copy.getEdges()) {
			edgeMap.put(heds.getEdge(e.getIndex()), e);
		}
		for (F f : copy.getFaces()) {
			faceMap.put(heds.getFace(f.getIndex()), f);
		}
		return copy;
	}
	
	
	/**
	 * calculates the cone peek of the circle left of medialEdge
	 * @param p the point to store the coordinate
	 * @param medialEdge the medial edge on which's left side p is
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double[] calculateConePeek(E medialEdge, Map<F, Double> rhoMap, AdapterSet aSet) {
		if (medialEdge.getLeftFace() != null){
			F medialFace = medialEdge.getLeftFace();
			double radius = Math.exp(rhoMap.get(medialFace));
			double rq = radius * radius;
			double[] c = aSet.getD(Position2d.class, medialFace); 
			double[] cq = {c[0]*c[0], c[1]*c[1]};
			double x = 2 * c[0];
			double y = 2 * c[1];
			double z = -1 + cq[0] + cq[1] - rq;
			double w =  1 + cq[0] + cq[1] - rq;
			return new double[] {x, z, y, w};
		} else {
			double[] p2 = aSet.getD(Position2d.class, medialEdge.getStartVertex());
			double[] p1 = aSet.getD(Position2d.class, medialEdge.getTargetVertex());
			double a = p2[1] - p1[1];
			double b = p1[0]*p2[1] - p1[1]*p2[0];
			double c = p1[0] - p2[0];
			double d = p1[1]*p2[0] - p1[0]*p2[1];
			return new double[] {a, b, c, -d};
		}
	}
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void inverseStereographicProjection(HDS hds, double scale, AdapterSet a){
		for (V v : hds.getVertices()){
			double[] vxy = a.getD(Position2d.class, v);
			double x = vxy[0] / scale;
			double y = vxy[1] / scale;
			double nx = 2 * x;
			double ny = x*x + y*y - 1;
			double nz = 2 * y;
			double nw = ny + 2;
			double[] p = {nx, ny, nz, nw};
			a.set(Position.class, v, p);
		}
	}
	
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double[] baryCenter(
		HalfEdgeDataStructure<V, E, F> graph,
		AdapterSet a
	) {
		List<V> vertices = graph.getVertices();
		double[] result = {0.0, 0.0};
		for (V v : vertices){
			double[] p = a.getD(Position2d.class, v);
			result[0] += p[0];
			result[1] += p[1];
		}
		Rn.times(result, 1.0 / graph.numVertices(), result);
		return result;
	}
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double meanLength(
		HalfEdgeDataStructure<V, E, F> graph,
		AdapterSet a
	) {
		List<V> vertices = graph.getVertices();
		double result = 0;
		for (V v : vertices){
			double[] p = a.getD(Position2d.class, v);
			result += Rn.euclideanNorm(p);
		}
		return result / graph.numVertices();
	}

	
	private static  <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void normalizeBeforeProjection(
		HalfEdgeDataStructure<V, E, F> medial, 
		Map<F, Double> rhoMap,
		AdapterSet a
	){
		double[] offset = baryCenter(medial, a);
		for (F f : medial.getFaces()){
			double[] fxy = a.getD(Position2d.class, f);
			Rn.subtract(fxy, fxy, offset);
			a.set(Position.class, f, fxy);
		}
		for (V v : medial.getVertices()){
			double[] vxy = a.getD(Position2d.class, v);
			Rn.subtract(vxy, vxy, offset);
			a.set(Position.class, v, vxy);
		}
		double scale = meanLength(medial, a);
		double logScale = Math.log(scale);
		for (F f : medial.getFaces()){
			double rho = rhoMap.get(f) - logScale;
			rhoMap.put(f, rho);
			double[] fxy = a.getD(Position2d.class, f);
			Rn.times(fxy, 1 / scale, fxy);
			a.set(Position.class, f, fxy);
		}
		for (V v : medial.getVertices()){
			double[] vxy = a.getD(Position2d.class, v);
			Rn.times(vxy, 1 / scale, vxy);
			a.set(Position.class, v, vxy);
		}

	}
	
}
