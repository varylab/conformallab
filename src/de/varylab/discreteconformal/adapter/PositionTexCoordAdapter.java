package de.varylab.discreteconformal.adapter;

import geom3d.Point;
import de.jtem.halfedge.jreality.adapter.CoordinateAdapter2Ifs;
import de.varylab.discreteconformal.heds.CoVertex;

public class PositionTexCoordAdapter implements CoordinateAdapter2Ifs<CoVertex> {

	private HyperbolicModel
		model = HyperbolicModel.Klein;
	private boolean
		projective = true;
	
	public PositionTexCoordAdapter(boolean projective) {
		this.projective = projective;
	}

	public PositionTexCoordAdapter(HyperbolicModel model, boolean projective) {
		this.model = model;
		this.projective = projective;
	}
	
	
	public double[] getCoordinate(CoVertex v) {
		Point t = v.getTextureCoord();
		if (projective) {
			switch (model) {
				case Klein:
					return new double[] {t.x(), t.y(), 0.0, t.z()};
				case Poincaré: 
				default:
					return new double[] {t.x(), t.y(), 0.0, t.z() + 1};
				case Halfplane:
					return new double[] {t.y(), 1, 0.0, t.z() - t.x()};
			}
		} else {
			switch (model) {
				case Klein:
					return new double[] {t.x() / t.z(), t.y() / t.z(), 0, 1};
				case Poincaré: 
				default:
					return new double[] {t.x() / t.z() + 1, t.y() / t.z() + 1, 0, 1};
				case Halfplane:
					return new double[] {t.y() / t.z() - t.x(), 1 / t.z() - t.x(), 0, 1};
			}
		}
	}

	public AdapterType getAdapterType() {
		return AdapterType.VERTEX_ADAPTER;
	}

	
	public void setModel(HyperbolicModel model) {
		this.model = model;
	}
	
	public void setProjective(boolean projective) {
		this.projective = projective;
	}
	
}
