package de.varylab.discreteconformal.heds.adapter;

import geom3d.Point;
import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

@Position
public class TexCoordPositionAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

	private HyperbolicModel
		model = HyperbolicModel.Klein;
	private boolean
		projective = true;
	
	public TexCoordPositionAdapter(boolean projective) {
		super(CoVertex.class, null, null, double[].class, true, false);
		this.projective = projective;
	}

	public TexCoordPositionAdapter(HyperbolicModel model, boolean projective) {
		this(projective);
		this.model = model;
	}
	
	
	@Override
	public double[] getVertexValue(CoVertex v, AdapterSet a) {
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
					return new double[] {t.x() / (t.z() + 1), t.y() / (t.z() + 1), 0, 1};
				case Halfplane:
					return new double[] {t.y() / (t.z() - t.x()), 1 / (t.z() - t.x()), 0, 1};
			}
		}
	}

	@Override
	public double getPriority() {
		return 10;
	}
	
	public void setModel(HyperbolicModel model) {
		this.model = model;
	}
	
	public void setProjective(boolean projective) {
		this.projective = projective;
	}
	
}
