package de.varylab.discreteconformal.heds.adapter;

import static de.jtem.halfedge.jreality.adapter.Adapter.AdapterType.VERTEX_ADAPTER;
import de.jtem.halfedge.jreality.adapter.LabelAdapter2Ifs;
import de.varylab.discreteconformal.heds.CVertex;

public class VertexLabelAdapter implements LabelAdapter2Ifs<CVertex> {

	public String getLabel(CVertex v) {
		return v.getIndex() + "";
	}

	public AdapterType getAdapterType() {
		return VERTEX_ADAPTER;
	}

}
