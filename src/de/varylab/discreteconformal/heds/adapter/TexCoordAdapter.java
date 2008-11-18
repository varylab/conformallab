package de.varylab.discreteconformal.heds.adapter;

import geom3d.Point;
import de.jtem.halfedge.jreality.adapter.TextCoordsAdapter2Ifs;
import de.varylab.discreteconformal.heds.CVertex;

public class TexCoordAdapter implements TextCoordsAdapter2Ifs<CVertex> {

	private boolean
		useProjectiveMap = true;
	
	
	public TexCoordAdapter() {

	}
	
	public TexCoordAdapter(boolean useProjective) {
		this.useProjectiveMap = useProjective;
	}
	
	
	public double[] getTextCoordinate(CVertex v) {
		Point t = v.getTextureCoord();
		if (useProjectiveMap) {
			return new double[] {t.x(), t.y(), 0.0, t.z()};
		} else {
			return new double[] {t.x() / t.z(), t.y() / t.z()};
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
