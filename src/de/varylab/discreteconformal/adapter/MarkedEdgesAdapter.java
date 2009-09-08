package de.varylab.discreteconformal.adapter;

import static de.jtem.halfedge.jreality.adapter.Adapter.AdapterType.EDGE_ADAPTER;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.jreality.adapter.ColorAdapter2Ifs;
import de.jtem.halfedge.jreality.adapter.RelRadiusAdapter2Ifs;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class MarkedEdgesAdapter
<V extends Vertex<V,E,F>,
E extends Edge<V,E,F>,
F extends Face<V,E,F>>  implements ColorAdapter2Ifs<Edge<?,?,?>>, RelRadiusAdapter2Ifs<Edge<?,?,?>> {

	private CuttingInfo<V, E, F>
		context = new CuttingInfo<V, E, F>();
	private Random
		rnd = new Random();
	private double[]
	    normalColor = {1, 1, 1};
	private Map<Set<E>, double[]>
		pathColors = new HashMap<Set<E>, double[]>();
	
	
	public MarkedEdgesAdapter() {
	}
	
	public MarkedEdgesAdapter(CuttingInfo<V, E, F> context) {
		this.context = context;
		updatePathColors();
	}
	
	public void setContext(CuttingInfo<V, E, F> context) {
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
		for (Set<E> path : context.paths) {
			rnd.setSeed(i);
			double[] color = new double[] {rnd.nextDouble(), rnd.nextDouble(), rnd.nextDouble()};
			pathColors.put(path, color);
			i++;
		}
	}
	
	
	@Override
	public double[] getColor(Edge<?,?,?> e) {
		for (Set<E> path : context.paths) {
			Set<E> coPath = context.pathCutMap.get(path);
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
	public double getReelRadius(Edge<?,?,?> e) {
		for (Set<E> path : context.paths) {
			Set<E> coPath = context.pathCutMap.get(path);
			if (path.contains(e) || path.contains(e.getOppositeEdge())) {
				return 3.0;
			}
			if(coPath != null) {
				if (coPath.contains(e) || coPath.contains(e.getOppositeEdge())) {
					return 3.0;
				}
			}
		}
		return 0.5;
	}
	
	
	@Override
	public AdapterType getAdapterType() {
		return EDGE_ADAPTER;
	}

}
