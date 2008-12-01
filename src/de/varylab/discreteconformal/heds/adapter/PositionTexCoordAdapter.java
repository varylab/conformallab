package de.varylab.discreteconformal.heds.adapter;

import geom3d.Point;
import de.jtem.halfedge.jreality.adapter.CoordinateAdapter2Ifs;
import de.varylab.discreteconformal.heds.CoVertex;

public class PositionTexCoordAdapter implements CoordinateAdapter2Ifs<CoVertex> {

	public double[] getCoordinate(CoVertex v) {
		Point t = v.getTextureCoord();
		return new double[] {t.x() / t.z(), t.y() / t.z(), 0.0};
	}

	public AdapterType getAdapterType() {
		return AdapterType.VERTEX_ADAPTER;
	}

}
