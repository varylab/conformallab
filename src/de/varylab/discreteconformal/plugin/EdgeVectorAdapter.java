package de.varylab.discreteconformal.plugin;

import java.util.List;
import java.util.Set;

import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.VectorField;
import de.jtem.halfedgetools.adapter.type.generic.EdgeVector;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

@VectorField
public class EdgeVectorAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

	private List<CoEdge>
		selectedEdged = null;
	private String
		name = "Edge Vectors";
	
	public EdgeVectorAdapter(List<CoEdge> path, String name) {
		super(null, CoEdge.class, null, double[].class, true, false);
		this.selectedEdged = path;
		this.name = name;
	}
	
	@Override
	public double[] getEdgeValue(CoEdge e, AdapterSet a) {
		if (selectedEdged.contains(e)) {
			return a.getD(EdgeVector.class, e);
		}
		return null;
	}
	
	@Override
	public String toString() {
		return name;
	}
	
}
