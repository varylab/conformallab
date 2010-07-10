package de.varylab.discreteconformal.heds.adapter;

import java.util.Set;

import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Radius;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

@Radius
public class MarkedEdgesRadiusAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, Double> {

	private CuttingInfo<CoVertex, CoEdge, CoFace>
		context = new CuttingInfo<CoVertex, CoEdge, CoFace>();
	
	
	public MarkedEdgesRadiusAdapter() {
		super(null, CoEdge.class, null, Double.class, true, false);
	}
	
	public MarkedEdgesRadiusAdapter(CuttingInfo<CoVertex, CoEdge, CoFace> context) {
		this();
		this.context = context;
	}
	
	public void setContext(CuttingInfo<CoVertex, CoEdge, CoFace> context) {
		this.context = context;
	}	
	
	
	@Override
	public Double getEdgeValue(CoEdge e, AdapterSet a) {
		for (Set<CoEdge> path : context.paths.keySet()) {
			Set<CoEdge> coPath = context.pathCutMap.get(path);
			if (path.contains(e) || path.contains(e.getOppositeEdge())) {
				return 4.0;
			}
			if(coPath != null) {
				if (coPath.contains(e) || coPath.contains(e.getOppositeEdge())) {
					return 4.0;
				}
			}
		}
		return 0.0;
	}

	@Override
	public double getPriority() {
		return 1;
	}
	
}
