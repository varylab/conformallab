package de.varylab.discreteconformal.heds.adapter;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Color;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

@Color
public class MarkedEdgesColorAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

	private CuttingInfo<CoVertex, CoEdge, CoFace>
		context = new CuttingInfo<CoVertex, CoEdge, CoFace>();
	private Random
		rnd = new Random();
	private double[]
	    normalColor = {1, 1, 1};
	private Map<Set<CoEdge>, double[]>
		pathColors = new HashMap<Set<CoEdge>, double[]>();
	
	
	public MarkedEdgesColorAdapter() {
		super(null, CoEdge.class, null, double[].class, true, false);
	}
		
	public MarkedEdgesColorAdapter(CuttingInfo<CoVertex, CoEdge, CoFace> context) {
		this();
		this.context = context;
		updatePathColors();
	}
	
	public void setContext(CuttingInfo<CoVertex, CoEdge, CoFace> context) {
		this.context = context;
		updatePathColors();
	}	
	
	
	private void updatePathColors() {
		pathColors.clear();
		if(context == null){
			System.err.println("context not set in MarkedEdgesAdapter");
		} else if (context.paths == null) {
			System.err.println("invalid context in MarkedEdgesAdapter");
		}
		
		int i = 0;
		for (Set<CoEdge> path : context.paths) {
			rnd.setSeed(i);
			double[] color = new double[] {rnd.nextDouble(), rnd.nextDouble(), rnd.nextDouble()};
			pathColors.put(path, color);
			i++;
		}
	}
	
	
	@Override
	public double[] getEdgeValue(CoEdge e, AdapterSet a) {
		for (Set<CoEdge> path : context.paths) {
			Set<CoEdge> coPath = context.pathCutMap.get(path);
			if (path.contains(e) || path.contains(e.getOppositeEdge())) {
				return pathColors.get(path);
			}
			if(coPath != null) {
				if (coPath.contains(e) || coPath.contains(e.getOppositeEdge())) {
					return pathColors.get(path);
				}
			}
		}
		return normalColor;
	}
	
	@Override
	public double getPriority() {
		return 1;
	}
	
	
}
