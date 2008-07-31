package de.varylab.discreteconformal.heds.adapter;

import static de.jtem.halfedge.jReality.interfaces.Adapter.AdapterType.VERTEX_ADAPTER;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.jReality.interfaces.LabelAdapter2Ifs;

public class VertexLabelAdapter implements LabelAdapter2Ifs {

	@Override
	public String getLabel(Node<?, ?, ?> node) {
		return node.getIndex() + "";
	}

	@Override
	public AdapterType getAdapterType() {
		return VERTEX_ADAPTER;
	}

}
