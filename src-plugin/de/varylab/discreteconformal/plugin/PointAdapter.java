package de.varylab.discreteconformal.plugin;

import static de.jtem.halfedge.jreality.adapter.Adapter.AdapterType.VERTEX_ADAPTER;

import java.util.HashSet;
import java.util.Set;

import de.jtem.halfedge.jreality.adapter.ColorAdapter2Ifs;
import de.jtem.halfedge.jreality.adapter.RelRadiusAdapter2Ifs;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.util.CuttingUtility.CuttingInfo;

public class PointAdapter implements ColorAdapter2Ifs<CoVertex>, RelRadiusAdapter2Ifs<CoVertex>  {

	private Set<CoVertex>
		rootCopySet = new HashSet<CoVertex>();
	private double[]
	    colorNormal = {1.0, 1.0, 1.0},
		colorMarked = {1.0, 0.0, 0.0};
	
	public PointAdapter() {

	}
	
	public void setContext(CuttingInfo<CoVertex, CoEdge, CoFace> context) {
		rootCopySet = context.getCopies(context.cutRoot);
		rootCopySet.add(context.cutRoot);
	}
	

	@Override
	public double getReelRadius(CoVertex node) {
		if (rootCopySet.contains(node)) {
			return 1.0;
		} else {
			return 0;
		}
	}
	
	
	@Override
	public double[] getColor(CoVertex node) {
		if (rootCopySet.contains(node)) {
			return colorMarked;
		} else {
			return colorNormal;
		}
	}

	@Override
	public AdapterType getAdapterType() {
		return VERTEX_ADAPTER;
	}

}
