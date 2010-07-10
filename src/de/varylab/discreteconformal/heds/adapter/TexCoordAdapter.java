package de.varylab.discreteconformal.heds.adapter;

import geom3d.Point;
import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.TexCoordinate;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

@TexCoordinate
public class TexCoordAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

	private HyperbolicModel
		model = HyperbolicModel.Klein;
	private boolean
		projective = true;
	private int
		priority = 1;
	
	
	public TexCoordAdapter(int priority) {
		super(CoVertex.class, null, null, double[].class, true, true);
		this.projective = false;
		this.priority = priority;
	}
	
	public TexCoordAdapter(boolean projective) {
		super(CoVertex.class, null, null, double[].class, true, true);
		this.projective = projective;
	}
	
	public TexCoordAdapter(HyperbolicModel model, boolean projective) {
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
					return new double[] {t.x() / t.z(), t.y() / t.z()};
				case Poincaré: 
				default:
					return new double[] {t.x() / (t.z() + 1), t.y() / (t.z() + 1)};
				case Halfplane:
					return new double[] {t.y() / (t.z() - t.x()), 1 / (t.z() - t.x())};
			}
		}
	}
	
	
	@Override
	public void setVertexValue(CoVertex v, double[] value, AdapterSet a) {
		if (value.length == 2) {
			v.getTextureCoord().get()[0] = value[0];
			v.getTextureCoord().get()[1] = value[1];
			v.getTextureCoord().get()[2] = 1.0;
		} else {
			v.getTextureCoord().set(value);
		}
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
