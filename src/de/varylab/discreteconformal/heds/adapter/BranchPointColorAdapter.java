package de.varylab.discreteconformal.heds.adapter;

import java.util.HashSet;
import java.util.Set;

import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Color;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

@Color
public class BranchPointColorAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

	private Set<CoVertex>
		rootCopySet = new HashSet<CoVertex>();
	private Set<CoVertex>
		branchSet = new HashSet<CoVertex>();
	private double[]
	    colorNormal = {0.5, 0.5, 0.5},
	    colorMarked = {1.0, 0.0, 0.0},
		colorMarked2 = {0.0, 1.0, 0.0};
	
	public BranchPointColorAdapter() {
		super(CoVertex.class, null, null, double[].class, true, false);
	}
	
	public void setContext(CuttingInfo<CoVertex, CoEdge, CoFace> context) {
		rootCopySet = context.getCopies(context.cutRoot);
		branchSet = context.getBranchSet();
	}
	
	@Override
	public double[] getVertexValue(CoVertex v, AdapterSet a) {
		if (rootCopySet.contains(v)) {
			return colorMarked;
		} else if (branchSet.contains(v)) {
			return colorMarked2;
		} else {
			return colorNormal;
		}
	}

	@Override
	public double getPriority() {
		return 1;
	}
	
}
