package de.varylab.discreteconformal.adapter;

import static de.jtem.halfedge.jreality.adapter.Adapter.AdapterType.EDGE_ADAPTER;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import de.jtem.halfedge.jreality.adapter.ColorAdapter2Ifs;
import de.jtem.halfedge.jreality.adapter.RelRadiusAdapter2Ifs;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class MarkedEdgesAdapter implements ColorAdapter2Ifs<CoEdge>, RelRadiusAdapter2Ifs<CoEdge> {

	private CuttingInfo<CoVertex, CoEdge, CoFace>
		context = new CuttingInfo<CoVertex, CoEdge, CoFace>();
	private Random
		rnd = new Random();
	private double[]
	    normalColor = {1, 1, 1};
	private Map<Set<CoEdge>, double[]>
		pathColors = new HashMap<Set<CoEdge>, double[]>();
	
	
	public MarkedEdgesAdapter() {
	}
	
	public MarkedEdgesAdapter(CuttingInfo<CoVertex, CoEdge, CoFace> context) {
		this.context = context;
		updatePathColors();
	}
	
	public void setContext(CuttingInfo<CoVertex, CoEdge, CoFace> context) {
		this.context = context;
		updatePathColors();
	}	
	
	
	private void updatePathColors() {
		pathColors.clear();
		for (Set<CoEdge> path : context.paths) {
			rnd.setSeed(path.size());
			double[] color = new double[] {rnd.nextDouble(), rnd.nextDouble(), rnd.nextDouble()};
			pathColors.put(path, color);
		}
	}
	
	
	@Override
	public double[] getColor(CoEdge e) {
		for (Set<CoEdge> path : context.paths) {
			Set<CoEdge> coPath = context.pathCutMap.get(path);
			if (path.contains(e) || path.contains(e.getOppositeEdge())) {
				return pathColors.get(path);
			}
			if (coPath.contains(e) || coPath.contains(e.getOppositeEdge())) {
				return pathColors.get(path);
			}
		}
		return normalColor;
	}
	
	@Override
	public double getReelRadius(CoEdge e) {
		for (Set<CoEdge> path : context.paths) {
			Set<CoEdge> coPath = context.pathCutMap.get(path);
			if (path.contains(e) || path.contains(e.getOppositeEdge())) {
				return 1.0;
			}
			if (coPath.contains(e) || coPath.contains(e.getOppositeEdge())) {
				return 1.0;
			}
		}
		return 0.0;
	}
	
	
	@Override
	public AdapterType getAdapterType() {
		return EDGE_ADAPTER;
	}

}
