package de.varylab.discreteconformal.plugin;

import java.util.LinkedList;
import java.util.List;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.varylab.discreteconformal.unwrapper.SphericalNormalizerPETSc;

public class SphericalNormalizationPlugin extends AlgorithmPlugin {

	public SphericalNormalizationPlugin() {
	}

	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.Geometry;
	}
	
	@Override
	public String getAlgorithmName() {
		return "Spherical Normalization";
	}

	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds, AdapterSet a, HalfedgeInterface hi) throws Exception {
		List<V> include = new LinkedList<V>(hi.getSelection().getVertices(hds));
		if (include.isEmpty()) {
			include = hds.getVertices();
		}
		SphericalNormalizerPETSc.normalize(hds, include, a, Position4d.class, Position.class);
		hi.update();
	}

}
