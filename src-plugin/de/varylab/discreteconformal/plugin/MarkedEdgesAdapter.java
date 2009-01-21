package de.varylab.discreteconformal.plugin;

import static de.jtem.halfedge.jreality.adapter.Adapter.AdapterType.EDGE_ADAPTER;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import de.jtem.halfedge.jreality.adapter.ColorAdapter2Ifs;
import de.jtem.halfedge.jreality.adapter.RelRadiusAdapter2Ifs;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.unwrapper.CHyperbolicLayout.HyperbolicLayoutContext;

public class MarkedEdgesAdapter implements ColorAdapter2Ifs<CoEdge>, RelRadiusAdapter2Ifs<CoEdge> {

	private HyperbolicLayoutContext
		context = new HyperbolicLayoutContext();
	private Random
		rnd = new Random();
	private double[]
	    normalColor = {1, 1, 1};
	private Map<Set<CoEdge>, double[]>
		pathColors = new HashMap<Set<CoEdge>, double[]>();
	
	
	public MarkedEdgesAdapter() {
	}
	
	public MarkedEdgesAdapter(HyperbolicLayoutContext context) {
		this.context = context;
		updatePathColors();
	}
	
	public void setContext(HyperbolicLayoutContext context) {
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
		if (!HalfEdgeUtils.isInteriorEdge(e)) {
			return new double[]{0, 0, 0};
		} else {
			return normalColor;
		}
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
		if (!HalfEdgeUtils.isInteriorEdge(e)) {
			return 0.3;
		} else {
			return 0.1;
		}
	}
	
	
	@Override
	public AdapterType getAdapterType() {
		return EDGE_ADAPTER;
	}

}
