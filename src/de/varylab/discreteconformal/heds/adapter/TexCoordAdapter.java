package de.varylab.discreteconformal.heds.adapter;

import de.jtem.halfedge.Node;
import de.jtem.halfedge.jReality.interfaces.TextCoordsAdapter2Ifs;
import de.varylab.discreteconformal.heds.CVertex;

public class TexCoordAdapter implements TextCoordsAdapter2Ifs {

	@Override
	public double[] getTextCoordinate(Node<?, ?, ?> node) {
		CVertex v = (CVertex)node;
		return v.getTextureCoord().get();
	}

	public AdapterType getAdapterType() {
		return AdapterType.VERTEX_ADAPTER;
	}

}
