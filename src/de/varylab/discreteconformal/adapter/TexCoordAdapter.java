package de.varylab.discreteconformal.adapter;

import geom3d.Point;
import de.jtem.halfedge.jreality.adapter.TextCoordsAdapter2Ifs;
import de.varylab.discreteconformal.heds.CoVertex;

public class TexCoordAdapter implements TextCoordsAdapter2Ifs<CoVertex> {

	private HyperbolicModel
		model = HyperbolicModel.Klein;
	private boolean
		projective = true;
	
	public TexCoordAdapter(boolean projective) {
		this.projective = projective;
	}
	
	public TexCoordAdapter(HyperbolicModel model, boolean projective) {
		this.model = model;
		this.projective = projective;
	}
	
	
	public double[] getTextCoordinate(CoVertex v) {
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
					return new double[] {t.x() / t.z(), t.y() / t.z()};
				case Poincaré: 
				default:
					return new double[] {t.x() / (t.z() + 1), t.y() / (t.z() + 1)};
				case Halfplane:
					return new double[] {t.y() / (t.z() - t.x()), 1 / (t.z() - t.x())};
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
