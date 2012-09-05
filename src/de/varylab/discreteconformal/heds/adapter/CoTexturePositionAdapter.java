package de.varylab.discreteconformal.heds.adapter;

import de.jreality.math.Rn;
import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.InterpolationMethod;

@TexturePosition
public class CoTexturePositionAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

	private HyperbolicModel
		model = HyperbolicModel.Klein;
	private InterpolationMethod
		interpolationMethod = InterpolationMethod.Incircle;
	private int
		priority = 1;
	
	public CoTexturePositionAdapter() {
		this(1);
	}
	
	public CoTexturePositionAdapter(int priority) {
		super(CoVertex.class, null, CoFace.class, double[].class, true, true);
		this.priority = priority;
	}
	
	
	public CoTexturePositionAdapter(HyperbolicModel model, InterpolationMethod interpolationMethod) {
		this(1);
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
						return new double[] {t[0] / t[3], t[1] / t[3]};
					case Poincaré: 
					default:
						return new double[] {t[0] / (t[3] + 1), t[1] / (t[3] + 1)};
					case Halfplane:
						return new double[] {t[1] / (t[3] - t[0]), 1 / (t[3] - t[0])};
				}
			case Circumcircle:
				t = Rn.times(null, t[3], t);
			default:
			case Incircle:	
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
		double[] t = v.T;
		t[0] = value[0];
		t[1] = value[1];
		t[2] = value.length > 2 ? value[2] : 0.0;
		t[3] = value.length > 3 ? value[3] : 1.0;
	}
	
	
	@Override
	public double[] getFaceValue(CoFace f, AdapterSet a) {
		return f.T;
	}
	
	
	@Override
	public void setFaceValue(CoFace f, double[] value, AdapterSet a) {
		double[] t = f.T;
		t[0] = value[0];
		t[1] = value[1];
		t[2] = value.length > 2 ? value[2] : 0.0;
		t[3] = value.length > 3 ? value[3] : 1.0;
	}
	
	
	public void setModel(HyperbolicModel model) {
		this.model = model;
	}
	
	public void setInterpolationMethod(InterpolationMethod interpolationMethod) {
		this.interpolationMethod = interpolationMethod;
	}
	
	@Override
	public double getPriority() {
		return priority;
	}
	public void setPriority(int priority) {
		this.priority = priority;
	}
	

}
