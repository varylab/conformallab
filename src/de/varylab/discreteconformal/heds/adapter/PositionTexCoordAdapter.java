package de.varylab.discreteconformal.heds.adapter;

import de.jtem.halfedge.Node;
import de.jtem.halfedge.jReality.interfaces.CoordinateAdapter2Ifs;
import de.varylab.discreteconformal.heds.CVertex;

public class PositionTexCoordAdapter implements CoordinateAdapter2Ifs {

	@Override
	public double[] getCoordinate(Node<?, ?, ?> node) {
		CVertex v = (CVertex)node;
		return v.getTextureCoord().get();
	}

	public AdapterType getAdapterType() {
		return AdapterType.VERTEX_ADAPTER;
	}

}
