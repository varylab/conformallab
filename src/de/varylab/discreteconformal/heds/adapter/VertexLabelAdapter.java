package de.varylab.discreteconformal.heds.adapter;

import static de.jtem.halfedge.jreality.adapter.Adapter.AdapterType.VERTEX_ADAPTER;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.jreality.adapter.LabelAdapter2Ifs;

public class VertexLabelAdapter implements LabelAdapter2Ifs {

	public String getLabel(Node<?, ?, ?> node) {
		return node.getIndex() + "";
	}

	public AdapterType getAdapterType() {
		return VERTEX_ADAPTER;
	}

}
