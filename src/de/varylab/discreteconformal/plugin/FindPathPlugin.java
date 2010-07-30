package de.varylab.discreteconformal.plugin;

import java.util.Iterator;
import java.util.List;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.CalculatorException;
import de.jtem.halfedgetools.adapter.CalculatorSet;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeSelection;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.varylab.discreteconformal.util.Search;

public class FindPathPlugin extends AlgorithmPlugin {

	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.Selection;
	}

	@Override
	public String getAlgorithmName() {
		return "Find Path";
	}

	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds, CalculatorSet c, HalfedgeInterface hcp) throws CalculatorException {
		if (hcp.getSelection().getVertices().size() != 2) {
			throw new RuntimeException("Select two vertices for path finding");
		}
		Iterator<V> vit = hcp.getSelection().getVertices(hds).iterator();
		V v1 = vit.next();
		V v2 = vit.next();
	    List<E> path = Search.bFS(v1, v2, true);
	    HalfedgeSelection sel = new HalfedgeSelection();
		for (E e : path) {
			sel.setSelected(e, true);
			sel.setSelected(e.getOppositeEdge(), true);
		}
		hcp.setSelection(sel);
	}

}
