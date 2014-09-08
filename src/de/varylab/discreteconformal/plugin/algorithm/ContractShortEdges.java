package de.varylab.discreteconformal.plugin.algorithm;

import java.awt.EventQueue;
import java.lang.reflect.InvocationTargetException;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JOptionPane;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.jtem.halfedgetools.selection.Selection;

public class ContractShortEdges extends AlgorithmPlugin {

	private Double
		maxLength = null;
	
	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.Editing;
	}

	@Override
	public String getAlgorithmName() {
		return "Contract Short Edges";
	}

	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds, AdapterSet a, HalfedgeInterface hcp) throws Exception {
		Selection s = hcp.getSelection();
		Runnable r = new Runnable() {
			@Override
			public void run() {
				String result = JOptionPane.showInputDialog(getOptionParent(), "Max", 1E-6);
				if (result == null) {
					maxLength = null;
				} else {
					maxLength = Double.parseDouble(result);
				}
			}
		};
		try {
			EventQueue.invokeAndWait(r);
		} catch (InvocationTargetException e) {
			throw (Exception)e.getCause();
		}
		if (maxLength == null) {
			return;
		}
		List<E> contract = new LinkedList<>();
		for (E e : hds.getPositiveEdges()) {
			Double l = a.get(Length.class, e, Double.class);
			if (l < maxLength) {
				contract.add(e);
			}
		}
		for (E e : contract) {
			if (!e.isValid()) continue;
			V v = TopologyAlgorithms.collapseEdge(e);
			TopologyAlgorithms.removeDigonsAt(v);
		}
		hcp.update();
		hcp.setSelection(s);
	}
	
}
