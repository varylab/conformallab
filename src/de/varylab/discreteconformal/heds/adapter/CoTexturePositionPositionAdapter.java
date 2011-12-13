package de.varylab.discreteconformal.heds.adapter;

import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

@Position
public class CoTexturePositionPositionAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

	private HyperbolicModel
		model = HyperbolicModel.Klein;
	private boolean
		projective = true;
	
	public CoTexturePositionPositionAdapter() {
		super(CoVertex.class, null, null, double[].class, true, true);
	}
	
	public CoTexturePositionPositionAdapter(HyperbolicModel model, boolean projective) {
		this();
		this.model = model;
		this.projective = projective;
	}
	
	@Override
	public double[] getVertexValue(CoVertex v, AdapterSet a) {
		double[] t = v.T;
		if (projective) {
			switch (model) {
				case Klein:
					return new double[] {t[0], t[1], t[2], t[3]};
				case Poincaré: 
				default:
					return new double[] {t[0], t[1], 0.0, t[3] + 1};
				case Halfplane:
					return new double[] {t[1], 1, 0.0, t[3] - t[0]};
			}
		} else {
			switch (model) {
				case Klein:
					return new double[] {t[0] / t[3], t[1] / t[3], 0.0};
				case Poincaré: 
				default:
					return new double[] {t[0] / (t[3] + 1), t[1] / (t[3] + 1), 0.0};
				case Halfplane:
					return new double[] {t[1] / (t[3] - t[0]), 1 / (t[3] - t[0]), 0.0};
			}
		}
	}
	
	@Override
	public void setVertexValue(CoVertex v, double[] value, AdapterSet a) {
		if (value.length >= 2) {
			v.T[0] = value[0];
			v.T[1] = value[1];
			v.T[2] = 0.0;
			v.T[3] = 1.0;
		}
		if (value.length >= 3) {
			v.T[2] = value[2];
		}
		if (value.length >= 4) {
			v.T[3] = value[3];
		}
	}

	@Override
	public double getPriority() {
		return 10;
	}
	
	public void setProjective(boolean projective) {
		this.projective = projective;
	}
	
	public void setModel(HyperbolicModel model) {
		this.model = model;
	}
	
}
