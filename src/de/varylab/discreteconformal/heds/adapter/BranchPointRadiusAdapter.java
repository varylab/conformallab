package de.varylab.discreteconformal.heds.adapter;

import java.util.HashSet;
import java.util.Set;

import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Radius;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

@Radius
public class BranchPointRadiusAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, Double> {

	private Set<CoVertex>
		rootCopySet = new HashSet<CoVertex>();
	private Set<CoVertex>
		branchSet = new HashSet<CoVertex>();
	
	public BranchPointRadiusAdapter() {
		super(CoVertex.class, null, null, Double.class, true, false);
	}
	
	public void setContext(CuttingInfo<CoVertex, CoEdge, CoFace> context) {
		rootCopySet = context.getCopies(context.cutRoot);
		branchSet = context.getBranchSet();
	}

	@Override
	public Double getVertexValue(CoVertex v, AdapterSet a) {
		if (rootCopySet.contains(v) || branchSet.contains(v)) {
			return 1.0;
		} else {
			return 0.0;
		}
	}
	
	
	@Override
	public double getPriority() {
		return 1;
	}
	
}
