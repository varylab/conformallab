package de.varylab.discreteconformal.adapter;

import geom3d.Point;
import de.jtem.halfedge.jreality.adapter.CoordinateAdapter2Ifs;
import de.varylab.discreteconformal.heds.CoVertex;

public class PositionTexCoordAdapter implements CoordinateAdapter2Ifs<CoVertex> {

	private HyperbolicModel
		model = HyperbolicModel.Klein;
	
	public PositionTexCoordAdapter() {
	}

	public PositionTexCoordAdapter(HyperbolicModel model) {
		this.model = model;
	}
	
	
	public double[] getCoordinate(CoVertex v) {
		Point t = v.getTextureCoord();
		switch (model) {
			case Klein:
				return new double[] {t.x(), t.y(), 0.0, t.z()};
			case Poincaré: 
			default:
				return new double[] {t.x(), t.y(), 0.0, t.z() + 1};
			case Halfplane:
				return new double[] {t.y(), 1, 0.0, t.z() - t.x()};
		}
	}

	public AdapterType getAdapterType() {
		return AdapterType.VERTEX_ADAPTER;
	}

	
	public void setModel(HyperbolicModel model) {
		this.model = model;
	}
	
}
