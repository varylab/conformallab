package de.varylab.discreteconformal.plugin;

import java.util.Set;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.CalculatorException;
import de.jtem.halfedgetools.adapter.CalculatorSet;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.Search;

public class CutToDiskPlugin extends AlgorithmPlugin {

	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.Topology;
	}

	@Override
	public String getAlgorithmName() {
		return "Cut To Disk";
	}

	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds, CalculatorSet c, HalfedgeInterface hcp) throws CalculatorException {
		V rootVertex = hds.getVertex(0);
		if (!hcp.getSelection().getVertices().isEmpty()) {
			Set<V> vertices = hcp.getSelection().getVertices(hds);
			rootVertex = vertices.iterator().next();
		}
		CuttingUtility.cutManifoldToDisk(hds, rootVertex, new Search.DefaultWeightAdapter<E>());
		hcp.update();
	}

}
