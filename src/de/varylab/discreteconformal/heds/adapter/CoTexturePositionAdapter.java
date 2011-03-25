package de.varylab.discreteconformal.heds.adapter;

import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

@TexturePosition
public class CoTexturePositionAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

	private HyperbolicModel
		model = HyperbolicModel.Klein;
	private boolean
		projective = true;
	private int
		priority = 1;
	
	public CoTexturePositionAdapter() {
		this(1);
	}
	
	public CoTexturePositionAdapter(int priority) {
		super(CoVertex.class, null, null, double[].class, true, true);
		this.projective = false;
		this.priority = priority;
	}
	
	public CoTexturePositionAdapter(boolean projective) {
		this(1);
		this.projective = projective;
	}
	
	public CoTexturePositionAdapter(HyperbolicModel model, boolean projective) {
		this(projective);
		this.model = model;
	}
	
	
	@Override
	public double[] getVertexValue(CoVertex v, AdapterSet a) {
		double[] t = v.T;
		if (projective) {
			switch (model) {
				case Klein:
					return t;
				case Poincaré: 
				default:
					return new double[] {t[0], t[1], 0.0, t[3] + 1};
				case Halfplane:
					return new double[] {t[1], 1, 0.0, t[3] - t[0]};
			}
		} else {
			switch (model) {
				case Klein:
					return new double[] {t[0] / t[3], t[1] / t[3]};
				case Poincaré: 
				default:
					return new double[] {t[0] / (t[3] + 1), t[1] / (t[3] + 1)};
				case Halfplane:
					return new double[] {t[1] / (t[3] - t[0]), 1 / (t[3] - t[0])};
			}
		}
	}
	
	
	@Override
	public void setVertexValue(CoVertex v, double[] value, AdapterSet a) {
		double[] t = v.T;
		t[0] = value[0];
		t[1] = value[1];
		t[2] = value.length > 2 ? value[2] : 0.0;
		t[3] = value.length > 3 ? value[3] : 0.0;
	}
	
	
	public void setModel(HyperbolicModel model) {
		this.model = model;
	}
	
	public void setProjective(boolean projective) {
		this.projective = projective;
	}
	
	@Override
	public double getPriority() {
		return priority;
	}
	public void setPriority(int priority) {
		this.priority = priority;
	}
	

}
