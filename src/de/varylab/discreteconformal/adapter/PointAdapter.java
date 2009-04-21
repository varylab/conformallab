package de.varylab.discreteconformal.adapter;

import static de.jtem.halfedge.jreality.adapter.Adapter.AdapterType.VERTEX_ADAPTER;

import java.util.HashSet;
import java.util.Set;

import de.jtem.halfedge.jreality.adapter.ColorAdapter2Ifs;
import de.jtem.halfedge.jreality.adapter.RelRadiusAdapter2Ifs;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class PointAdapter implements ColorAdapter2Ifs<CoVertex>, RelRadiusAdapter2Ifs<CoVertex>  {

	private Set<CoVertex>
		rootCopySet = new HashSet<CoVertex>();
	private Set<CoVertex>
		branchSet = new HashSet<CoVertex>();
	private double[]
	    colorNormal = {1.0, 1.0, 1.0},
	    colorMarked = {1.0, 0.0, 0.0},
		colorMarked2 = {0.0, 1.0, 0.0};
	
	public PointAdapter() {

	}
	
	public void setContext(CuttingInfo<CoVertex, CoEdge, CoFace> context) {
		rootCopySet = context.getCopies(context.cutRoot);
		branchSet = context.getBranchSet();
	}
	

	@Override
	public double getReelRadius(CoVertex node) {
		if (rootCopySet.contains(node) || branchSet.contains(node)) {
			return 1.0;
		} else {
			return 0;
		}
	}
	
	
	@Override
	public double[] getColor(CoVertex node) {
		if (rootCopySet.contains(node)) {
			return colorMarked;
		} else if (branchSet.contains(node)) {
			return colorMarked2;
		} else {
			return colorNormal;
		}
	}

	@Override
	public AdapterType getAdapterType() {
		return VERTEX_ADAPTER;
	}

}
