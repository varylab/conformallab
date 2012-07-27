package de.varylab.discreteconformal.plugin;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.varylab.discreteconformal.unwrapper.koebe.KoebePolyhedron;
import de.varylab.discreteconformal.unwrapper.koebe.KoebePolyhedronContext;

public class KoebePolyhedronPlugin extends AlgorithmPlugin {

	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.Geometry;
	}

	@Override
	public String getAlgorithmName() {
		return "Koebe Polyhedron";
	}

	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds, AdapterSet a, HalfedgeInterface hi) {
		KoebePolyhedronContext<V, E, F> context = KoebePolyhedron.contructKoebePolyhedron(hds, a);
		hi.set(context.polyhedron);
	}

}
