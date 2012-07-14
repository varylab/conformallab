package de.varylab.discreteconformal.heds.adapter;

import de.jreality.math.Rn;
import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.InterpolationMethod;

@Position
public class CoTexturePositionPositionAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

	private HyperbolicModel
		model = HyperbolicModel.Klein;
	private InterpolationMethod
		interpolationMethod = InterpolationMethod.Optimal;
	
	public CoTexturePositionPositionAdapter() {
		super(CoVertex.class, null, null, double[].class, true, true);
	}
	
	public CoTexturePositionPositionAdapter(HyperbolicModel model, InterpolationMethod interpolationMethod) {
		this();
		this.model = model;
		this.interpolationMethod = interpolationMethod;
	}
	
	@Override
	public double[] getVertexValue(CoVertex v, AdapterSet a) {
		double[] t = v.T;
		switch (interpolationMethod) {
			case Linear:
				switch (model) {
					case Klein:
						return new double[] {t[0] / t[3], t[1] / t[3], t[2] / t[3], 1};
					case Poincaré: 
					default:
						return new double[] {t[0] / (t[3] + 1), t[1] / (t[3] + 1), 0, 1};
					case Halfplane:
						return new double[] {t[1] / (t[3] - t[0]), 1 / (t[3] - t[0]), 0, 1};
				}
			case Circumcircle:
				t = Rn.times(null, t[3], t);
			default:
			case Optimal:	
				switch (model) {
					case Klein:
						return t;
					case Poincaré: 
					default:
						return new double[] {t[0], t[1], 0.0, t[3] + 1};
					case Halfplane:
						return new double[] {t[1], 1, 0.0, t[3] - t[0]};
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
	
	public void setInterpolationMethod(InterpolationMethod interpolationMethod) {
		this.interpolationMethod = interpolationMethod;
	}
	
	public void setModel(HyperbolicModel model) {
		this.model = model;
	}
	
}
