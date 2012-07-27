package de.varylab.discreteconformal.unwrapper.koebe;

import java.util.HashMap;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;

public class KoebePolyhedronContext<
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>
> {
	
	public HalfEdgeDataStructure<V, E, F>
		circlePattern = null,
		medial = null,
		polyhedron = null;
	public HashMap<E, E>
		edgeEdgeMap = null;
	public HashMap<E, V>
		edgeVertexMap = null;
	public HashMap<F, F>
		faceFaceMap = null;
	public HashMap<V, F>
		vertexFaceMap = new HashMap<V, F>();
	public V
		northPole = null;
	
	public HashMap<E, E> getEdgeEdgeMap() {
		return edgeEdgeMap;
	}
	public HalfEdgeDataStructure<V, E, F> getMedial() {
		return medial;
	}
	public HalfEdgeDataStructure<V, E, F> getPolyeder() {
		return polyhedron;
	}
	public HashMap<E, V> getEdgeVertexMap() {
		return edgeVertexMap;
	}
	public V getNorthPole() {
		return northPole;
	}
	public void setNorthPole(V northPole) {
		this.northPole = northPole;
	}
}