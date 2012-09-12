package de.varylab.discreteconformal.unwrapper.quasiisothermic;

import java.util.HashMap;
import java.util.Map;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;

public class DBFSolution <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>,
	HDS extends HalfEdgeDataStructure<V, E, F>
> {

	public HDS
		hds = null;
	public Map<E, Double>
		solutionAlphaMap = new HashMap<E, Double>(),
		initialAlphaMap = new HashMap<E, Double>(); 
	
}
