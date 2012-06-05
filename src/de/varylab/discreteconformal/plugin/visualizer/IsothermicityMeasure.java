package de.varylab.discreteconformal.plugin.visualizer;


import static java.lang.Math.sin;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.CurvatureFieldMin;
import de.jtem.halfedgetools.adapter.type.Normal;
import de.jtem.halfedgetools.adapter.type.generic.EdgeVector;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.Plugin;
import de.varylab.discreteconformal.unwrapper.isothermic.IsothermicUtility;

public class IsothermicityMeasure extends Plugin {


	@Retention(RetentionPolicy.RUNTIME)
	@Target(ElementType.TYPE)
	public static @interface CurvatureAngle {}
	
	@CurvatureAngle
	public class CurvatureAngleAdapter extends AbstractAdapter<Double> {

		public CurvatureAngleAdapter() {
			super(Double.class, true, false);
		}
		
		@Override
		public <
			V extends Vertex<V,E,F>, 
			E extends Edge<V,E,F>, 
			F extends Face<V,E,F>
		> Double getE(E e, AdapterSet a) {
			double[] N = a.getD(Normal.class, e);
			double[] Kmin = a.getD(CurvatureFieldMin.class, e);
			double[] E = a.getD(EdgeVector.class, e);
			return IsothermicUtility.getSignedAngle(N, Kmin, E);
		}
		
		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return Edge.class.isAssignableFrom(nodeClass);
		}
		
		@Override
		public String toString() {
			return "Principle Curvature Angle";
		}
		
	}
	
	
	public class SinConditionAdapter extends AbstractAdapter<Double> {
		
		public SinConditionAdapter() {
			super(Double.class, true, false);
		}
		
		@Override
		public <
			V extends Vertex<V, E, F>,
			E extends Edge<V, E, F>,
			F extends Face<V, E, F>
		> Double getV(V v, AdapterSet a) {
			double sl = 1;
			double sr = 1;
			if (HalfEdgeUtils.isBoundaryVertex(v)) {
				return 0.0;
			}
			for (E e : HalfEdgeUtils.incomingEdges(v)) {
				double alphaE = a.get(CurvatureAngle.class, e, Double.class);
				double alphaEPrev = a.get(CurvatureAngle.class, e.getPreviousEdge(), Double.class);
				double alphaENext = a.get(CurvatureAngle.class, e.getNextEdge(), Double.class);
				double alphaRight = IsothermicUtility.calculateTriangleAngle(alphaE, alphaEPrev, alphaENext);
				double alphaLeft = IsothermicUtility.calculateTriangleAngle(alphaENext, alphaEPrev, alphaE);
				sr *= sin(alphaRight);
				sl *= sin(alphaLeft);
			}
			return (sl - sr)*(sl - sr);
		}
		
		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return Vertex.class.isAssignableFrom(nodeClass);
		}
		
		@Override
		public String toString() {
			return "Sinus-Condition";
		}
		
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		c.getPlugin(HalfedgeInterface.class).addAdapter(new CurvatureAngleAdapter(), true);
		c.getPlugin(HalfedgeInterface.class).addAdapter(new SinConditionAdapter(), true);
	}
	
}
