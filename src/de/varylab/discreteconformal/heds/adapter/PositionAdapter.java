package de.varylab.discreteconformal.heds.adapter;

import static de.jtem.halfedge.jreality.adapter.Adapter.AdapterType.VERTEX_ADAPTER;
import geom3d.Point;
import de.jtem.halfedge.jreality.adapter.CoordinateAdapter2Heds;
import de.jtem.halfedge.jreality.adapter.CoordinateAdapter2Ifs;
import de.varylab.discreteconformal.heds.CoVertex;

public class PositionAdapter implements CoordinateAdapter2Heds<CoVertex>,
		CoordinateAdapter2Ifs<CoVertex> {

	public void setCoordinate(CoVertex v, double[] coord) {
		v.getPosition().set(coord);
	}

	public double[] getCoordinate(CoVertex v) {
		Point p = v.getPosition();
		System.out.println(p);
		return new double[] {p.x(), p.y(), p.z(), 1.0};
	}

	public AdapterType getAdapterType() {
		return VERTEX_ADAPTER;
	}
	
}
