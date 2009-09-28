package de.varylab.discreteconformal.adapter;

import static de.jtem.halfedgetools.jreality.adapter.Adapter.AdapterType.VERTEX_ADAPTER;
import de.jtem.halfedgetools.jreality.adapter.LabelAdapter2Ifs;
import de.varylab.discreteconformal.heds.CoVertex;

public class VertexLabelAdapter implements LabelAdapter2Ifs<CoVertex> {

	public String getLabel(CoVertex v) {
		return v.getIndex() + "";
	}

	public AdapterType getAdapterType() {
		return VERTEX_ADAPTER;
	}

}
