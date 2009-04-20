package de.varylab.discreteconformal.heds.adapter;

import geom3d.Point;
import de.jtem.halfedge.jreality.adapter.CoordinateAdapter2Ifs;
import de.varylab.discreteconformal.heds.CoVertex;

public class PositionTexCoordAdapter implements CoordinateAdapter2Ifs<CoVertex> {

	private boolean
		poincare = false;
	
	public PositionTexCoordAdapter() {
	}

	public PositionTexCoordAdapter(boolean poincare) {
		this.poincare = poincare;
	}
	
	
	public void setPoincare(boolean poincare) {
		this.poincare = poincare;
	}
	
	public double[] getCoordinate(CoVertex v) {
		Point t = v.getTextureCoord();
		if (poincare) {
			return new double[] {t.x(), t.y(), 0.0, t.z() + 1};
		} else {
			return new double[] {t.x(), t.y(), 0.0, t.z()};
		}
	}

	public AdapterType getAdapterType() {
		return AdapterType.VERTEX_ADAPTER;
	}

}
