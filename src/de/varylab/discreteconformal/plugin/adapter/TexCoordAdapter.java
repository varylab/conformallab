package de.varylab.discreteconformal.plugin.adapter;

import geom3d.Point;
import de.jtem.halfedge.jreality.adapter.TextCoordsAdapter2Ifs;
import de.varylab.discreteconformal.heds.CoVertex;

public class TexCoordAdapter implements TextCoordsAdapter2Ifs<CoVertex> {

	private boolean
		useProjectiveMap = true,
		poincare = false;
	
	
	public TexCoordAdapter() {

	}
	
	public TexCoordAdapter(boolean useProjective, boolean poincare) {
		this.useProjectiveMap = useProjective;
		this.poincare = poincare;
	}
	
	public TexCoordAdapter(boolean useProjective) {
		this.useProjectiveMap = useProjective;
	}
	
	public double[] getTextCoordinate(CoVertex v) {
		Point t = v.getTextureCoord();
		if (poincare) {
			return new double[] {t.x(), t.y(), 0.0, t.z() + 1.0};
		} else {
			if (useProjectiveMap) {
				return new double[] {t.x(), t.y(), 0.0, t.z()};
			} else {
				return new double[] {t.x(), t.y(), 0.0, t.z()};
			}
		}
	}

	public AdapterType getAdapterType() {
		return AdapterType.VERTEX_ADAPTER;
	}

	public boolean isUseProjectiveMap() {
		return useProjectiveMap;
	}

	public void setUseProjectiveMap(boolean useProjectiveMap) {
		this.useProjectiveMap = useProjectiveMap;
	}

}
