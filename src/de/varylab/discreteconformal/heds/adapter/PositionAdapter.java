package de.varylab.discreteconformal.heds.adapter;

import static de.jtem.halfedge.jreality.adapter.Adapter.AdapterType.VERTEX_ADAPTER;
import geom3d.Point;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.jreality.adapter.CoordinateAdapter2Heds;
import de.jtem.halfedge.jreality.adapter.CoordinateAdapter2Ifs;
import de.varylab.discreteconformal.heds.CVertex;

public class PositionAdapter implements CoordinateAdapter2Heds,
		CoordinateAdapter2Ifs {

	public void setCoordinate(Node<?, ?, ?> node, double[] coord) {
		CVertex v = (CVertex)node;
		v.getPosition().set(coord);
	}

	public double[] getCoordinate(Node<?, ?, ?> node) {
		CVertex v = (CVertex)node;
		Point p = v.getPosition();
		return new double[] {p.x(), p.y(), p.z(), 1.0};
	}

	public AdapterType getAdapterType() {
		return VERTEX_ADAPTER;
	}
	
}
