package de.varylab.discreteconformal.heds.adapter;

import static de.jtem.halfedge.jreality.adapter.Adapter.AdapterType.VERTEX_ADAPTER;
import geom3d.Point;
import de.jtem.halfedge.jreality.adapter.CoordinateAdapter2Heds;
import de.jtem.halfedge.jreality.adapter.CoordinateAdapter2Ifs;
import de.varylab.discreteconformal.heds.CVertex;

public class PositionAdapter implements CoordinateAdapter2Heds<CVertex>,
		CoordinateAdapter2Ifs<CVertex> {

	public void setCoordinate(CVertex v, double[] coord) {
		v.getPosition().set(coord);
	}

	public double[] getCoordinate(CVertex v) {
		Point p = v.getPosition();
		return new double[] {p.x(), p.y(), p.z(), 1.0};
	}

	public AdapterType getAdapterType() {
		return VERTEX_ADAPTER;
	}
	
}
