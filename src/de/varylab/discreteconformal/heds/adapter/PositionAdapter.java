package de.varylab.discreteconformal.heds.adapter;

import static de.jtem.halfedge.jReality.interfaces.Adapter.AdapterType.VERTEX_ADAPTER;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.jReality.interfaces.CoordinateAdapter2Heds;
import de.jtem.halfedge.jReality.interfaces.CoordinateAdapter2Ifs;
import de.varylab.discreteconformal.heds.CVertex;

public class PositionAdapter implements CoordinateAdapter2Heds,
		CoordinateAdapter2Ifs {

	public void setCoordinate(Node<?, ?, ?> node, double[] coord) {
		CVertex v = (CVertex)node;
		v.getPosition().set(coord);
	}

	public double[] getCoordinate(Node<?, ?, ?> node) {
		CVertex v = (CVertex)node;
		return v.getPosition().get();
	}

	public AdapterType getAdapterType() {
		return VERTEX_ADAPTER;
	}
	
}
