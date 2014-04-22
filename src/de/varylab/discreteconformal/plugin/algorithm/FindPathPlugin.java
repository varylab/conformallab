package de.varylab.discreteconformal.plugin.algorithm;

import java.util.Iterator;
import java.util.List;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.SelectionInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.jtem.halfedgetools.selection.Selection;
import de.jtem.jrworkspace.plugin.Controller;
import de.varylab.discreteconformal.util.Search;

public class FindPathPlugin extends AlgorithmPlugin {

	private final Integer
		PATH_CHANNEL = 293842934;
	
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
	> void execute(HDS hds, AdapterSet a, HalfedgeInterface hcp) {
		if (hcp.getSelection().getVertices().size() != 2) {
			throw new RuntimeException("Select two vertices for path finding");
		}
		Iterator<V> vit = hcp.getSelection().getVertices(hds).iterator();
		V v1 = vit.next();
		V v2 = vit.next();
	    List<E> path = Search.bFS(v1, v2, true);
	    Selection sel = new Selection();
		for (E e : path) {
			sel.add(e, PATH_CHANNEL);
			sel.add(e.getOppositeEdge(), PATH_CHANNEL);
		}
		hcp.addSelection(sel);
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		c.getPlugin(SelectionInterface.class).registerChannelName(PATH_CHANNEL, "Path Edges");
	}

}
