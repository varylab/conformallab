package de.varylab.discreteconformal.adapter;

import geom3d.Point;
import de.jtem.halfedge.jreality.adapter.TextCoordsAdapter2Ifs;
import de.varylab.discreteconformal.heds.CoVertex;

public class TexCoordAdapter implements TextCoordsAdapter2Ifs<CoVertex> {

	private HyperbolicModel
		model = HyperbolicModel.Klein;
	
	public TexCoordAdapter() {
	}
	
	public TexCoordAdapter(HyperbolicModel model) {
		this.model = model;
	}
	
	
	public double[] getTextCoordinate(CoVertex v) {
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
