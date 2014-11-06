package de.varylab.discreteconformal.plugin.visualizer;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.LengthTex;
import de.jtem.halfedgetools.plugin.data.DataSourceProvider;
import de.jtem.jrworkspace.plugin.Plugin;

public class DiscreteConformalEquivalencemMeasure extends Plugin implements DataSourceProvider {

	private class DataSource extends AbstractAdapter<Double> {

		public DataSource() {
			super(Double.class, true, false);
		}

		@Override
		public <
			V extends Vertex<V, E, F>, 
			E extends Edge<V, E, F>, 
			F extends Face<V, E, F>
		> Double getE(E e, AdapterSet a) {
//			Integer selectionChannel = a.get(Selection.class, e, Integer.class);
//			if (selectionChannel == null) {
//				return 0.0;
//			}
			if (HalfEdgeUtils.isBoundaryEdge(e)) {
				return 0.0;
			}
			double lcr = 1.0;
			double tlcr = 1.0;
			if (e.isPositive()) {
				e = e.getOppositeEdge();
			}
			if (e.getLeftFace() != null) {
				E next = e.getNextEdge();
				E prev = e.getPreviousEdge();
				double ln = a.get(Length.class, next, Double.class).doubleValue();
				double lp = a.get(Length.class, prev, Double.class).doubleValue();
				double ltn = a.get(LengthTex.class, next, Double.class).doubleValue();
				double ltp = a.get(LengthTex.class, prev, Double.class).doubleValue();
				lcr *= ln / lp;
				tlcr *= ltn / ltp;
			}
			E o = e.getOppositeEdge();
			if (o.getLeftFace() != null) {
				E next = o.getNextEdge();
				E prev = o.getPreviousEdge();
				double ln = a.get(Length.class, next, Double.class).doubleValue();
				double lp = a.get(Length.class, prev, Double.class).doubleValue();
				double ltn = a.get(LengthTex.class, next, Double.class).doubleValue();
				double ltp = a.get(LengthTex.class, prev, Double.class).doubleValue();
				lcr *= ln / lp;
				tlcr *= ltn / ltp;				
			}
			return Math.abs(lcr - tlcr);
		}
		
		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return Edge.class.isAssignableFrom(nodeClass);
		}

		@Override
		public String toString() {
			return "Conformal Equivalence";
		}
		
	}
	
	@Override
	public AdapterSet getDataSources() {
		return new AdapterSet(new DataSource());
	}

}
