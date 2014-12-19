package de.varylab.discreteconformal.plugin.algorithm;

import java.util.List;
import java.util.Set;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.jtem.halfedgetools.selection.Selection;
import de.varylab.discreteconformal.util.HomotopyUtility;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

public class SelectHomotopyGeneratorsPlugin extends AlgorithmPlugin {

	private final int
		CHANNEL_OFFSET = 347434567;
	
	@Override
	public String getAlgorithmName() {
		return "Homotopy Generators";
	}

	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.Selection;
	}
	
	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds, final AdapterSet a, HalfedgeInterface hi) throws Exception {
		Set<V> vSelection = hi.getSelection().getVertices(hds);
		if (vSelection.size() != 1) {
			throw new IllegalArgumentException("Select a unique root vertex");
		}
		V root = vSelection.iterator().next();
		WeightAdapter<E> lengthWeights = new WeightAdapter<E>() {
			@Override
			public double getWeight(E e) {
				return a.get(Length.class, e, Double.class).doubleValue();
			}
		};
		List<Set<E>> paths = HomotopyUtility.getGeneratorPaths(root, lengthWeights);
		Selection s = new Selection();
		int index = 0;
		for (Set<E> path : paths) {
			s.addAll(path, CHANNEL_OFFSET + index++);
		}
		hi.addSelection(s);
	}

}
