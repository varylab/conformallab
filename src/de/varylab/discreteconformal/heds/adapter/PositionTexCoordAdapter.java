package de.varylab.discreteconformal.heds.adapter;

import geom3d.Point;
import de.jtem.halfedge.jreality.adapter.CoordinateAdapter2Ifs;
import de.varylab.discreteconformal.heds.CVertex;

public class PositionTexCoordAdapter implements CoordinateAdapter2Ifs<CVertex> {

	public double[] getCoordinate(CVertex v) {
		Point t = v.getTextureCoord();
		return new double[] {t.x() / t.z(), t.y() / t.z(), 0.0};
	}

	public AdapterType getAdapterType() {
		return AdapterType.VERTEX_ADAPTER;
	}

}
