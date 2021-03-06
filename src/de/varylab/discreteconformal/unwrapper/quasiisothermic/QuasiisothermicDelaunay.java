package de.varylab.discreteconformal.unwrapper.quasiisothermic;

import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryEdge;
import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.util.ConsistencyCheck;
import de.jtem.halfedgetools.util.TriangulationException;
import de.varylab.discreteconformal.heds.adapter.types.OppositeAngle;


/**
 * Construct a delaunay triangulation from a given triangulation
 * <p>
 * Copyright 2005 <a href="http://www.sechel.de">Stefan Sechelmann</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Stefan Sechelmann
 */
public class QuasiisothermicDelaunay {

	
	/**
	 * Checks if this edge is locally delaunay
	 * @param edge the edge to check
	 * @return the check result
	 * @throws TriangulationException
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> boolean isDelaunay(E edge, AdapterSet a) throws TriangulationException{
		Double gamma = a.get(OppositeAngle.class, edge, Double.class);
		Double delta = a.get(OppositeAngle.class, edge.getOppositeEdge(), Double.class);
		return gamma + delta <= Math.PI;
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V,E,F>
	> boolean isDelaunay(HDS graph, AdapterSet a) {
		for(E edge : graph.getEdges()){
			try{
				if(!isDelaunay(edge, a))
				return false;
			} catch (TriangulationException e) {
				return false;
			}
		}
		return true;
	}
	
	
	/**
	 * Returns the positive edges belonging to the kite of the edge 
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<E> getPositiveKiteBorder(E edge){
		LinkedList<E> result = new LinkedList<E>();
		E e1 = edge.getNextEdge();
		E e2 = edge.getPreviousEdge();
		E e3 = edge.getOppositeEdge().getNextEdge();
		E e4 = edge.getOppositeEdge().getPreviousEdge();
		if (!e1.isPositive())
			e1 = e1.getOppositeEdge();
		if (!e2.isPositive())
			e2 = e2.getOppositeEdge();
		if (!e3.isPositive())
			e3 = e3.getOppositeEdge();
		if (!e4.isPositive())
			e4 = e4.getOppositeEdge();
		result.add(e1);
		result.add(e2);
		result.add(e3);
		result.add(e4);
		return result;
	}
	
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void flip(E e, AdapterSet a) throws TriangulationException{
		F leftFace = e.getLeftFace();
		F rightFace = e.getRightFace();
		if (leftFace == rightFace) {
			return;
		}
		// calculate new angles
//		System.out.println("flipped edge " + e);
		setMetricFlipAngles(e, a);
//		setPtolemyFlipAngles(e, a); // does not terminate
		
		E a1 = e.getOppositeEdge().getNextEdge();
		E a2 = a1.getNextEdge();
		E b1 = e.getNextEdge();
		E b2 = b1.getNextEdge();
		
		V v1 = e.getStartVertex();
		V v2 = e.getTargetVertex();
		V v3 = a1.getTargetVertex();
		V v4 = b1.getTargetVertex();

		//new connections
		e.linkNextEdge(a2);
		e.linkPreviousEdge(b1);
		e.getOppositeEdge().linkNextEdge(b2);
		e.getOppositeEdge().linkPreviousEdge(a1);
		e.setTargetVertex(v3);
		e.getOppositeEdge().setTargetVertex(v4);
		
		a2.linkNextEdge(b1);
		b2.linkNextEdge(a1);
		
		//set faces
		b2.setLeftFace(rightFace);
		a2.setLeftFace(leftFace);
		
		b2.setTargetVertex(v1);
		a2.setTargetVertex(v2);
		a1.setTargetVertex(v3);
		b1.setTargetVertex(v4);
	}
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void setMetricFlipAngles(E e, AdapterSet a) throws TriangulationException{
		F leftFace = e.getLeftFace();
		F rightFace = e.getRightFace();
		if (leftFace == rightFace) {
			return;
		}
		E a1 = e.getOppositeEdge().getNextEdge();
		E a2 = a1.getNextEdge();
		E b1 = e.getNextEdge();
		E b2 = b1.getNextEdge();
		
		// calculate new angles
		Double a1a = a.get(OppositeAngle.class, a1, Double.class); 
		Double a2a = a.get(OppositeAngle.class, a2, Double.class); 
		Double b1a = a.get(OppositeAngle.class, b1, Double.class); 
		Double b2a = a.get(OppositeAngle.class, b2, Double.class); 
		Double ea = a.get(OppositeAngle.class, e, Double.class);
		Double eoa = a.get(OppositeAngle.class, e.getOppositeEdge(), Double.class);
		
		double la1 = sin(ea)*sin(a1a)/(sin(b2a)*sin(eoa));
		double lep = sqrt(la1*la1 + 1 - 2*la1*cos(a2a + b1a));
		double g3 = sin(a2a + b1a) / lep;
		double g2 = PI - g3 - b1a - a2a;
		double g1 = ea - g2;
		double g4 = eoa - g3;
		
		a.set(OppositeAngle.class, e, b2a + a1a);
		a.set(OppositeAngle.class, e.getOppositeEdge(), b1a + a2a);
		a.set(OppositeAngle.class, a1, g2);
		a.set(OppositeAngle.class, a2, g1);
		a.set(OppositeAngle.class, b1, g4);
		a.set(OppositeAngle.class, b2, g3);
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void setPtolemyFlipAngles(E e, AdapterSet a) throws TriangulationException{
		F leftFace = e.getLeftFace();
		F rightFace = e.getRightFace();
		if (leftFace == rightFace) {
			return;
		}
		E eo = e.getOppositeEdge();
		E a1 = e.getOppositeEdge().getNextEdge();
		E a2 = a1.getNextEdge();
		E b1 = e.getNextEdge();
		E b2 = b1.getNextEdge();
		
		// calculate new angles
		Double a1a = a.get(OppositeAngle.class, a1, Double.class); 
		Double a2a = a.get(OppositeAngle.class, a2, Double.class); 
		Double b1a = a.get(OppositeAngle.class, b1, Double.class); 
		Double b2a = a.get(OppositeAngle.class, b2, Double.class); 
		Double ea = a.get(OppositeAngle.class, e, Double.class);
		Double eoa = a.get(OppositeAngle.class, e.getOppositeEdge(), Double.class);
		
		
		double le = sin(ea) / sin(b2a);
		double la1 = le * sin(a1a) / sin(eoa);
		double la2 = le * sin(a2a) / sin(eoa);
		double lb1 = le * sin(b1a) / sin(ea);
		double lb2 = 1.0;
		double lep = (la1*lb1 + la2*lb2) / le;
		
		a.set(OppositeAngle.class, e, cosineLaw(lep, lb1, la2));
		a.set(OppositeAngle.class, eo, cosineLaw(lep, lb2, la1));
		a.set(OppositeAngle.class, a1, cosineLaw(la1, lb2, lep));
		a.set(OppositeAngle.class, a2, cosineLaw(la2, lb1, lep));
		a.set(OppositeAngle.class, b1, cosineLaw(lb2, la1, lep));
		a.set(OppositeAngle.class, b2, cosineLaw(lb1, la2, lep));
	}
	
	
	/**
	 * Calculate the opposite angle to length a
	 * @param a
	 * @param b
	 * @param c
	 * @return
	 */
	protected static double cosineLaw(double a, double b, double c) {
		return acos((a*a - b*b - c*c) / (b*c));
	}
	
	
	/**
	 * Constructs the delaunay triangulation of the given structure.
	 * @param graph must be a triangulation
	 * @throws TriangulationException if the given graph is no triangulation or 
	 * if the trangle inequation doesn't hold for some triangle
	 * @return The new delaunay angles as adapter
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void constructDelaunay(HDS graph, AdapterSet a) throws TriangulationException {
		if (!ConsistencyCheck.isTriangulation(graph)) {
			throw new TriangulationException("Graph is no triangulation!");
		}
		HashSet<E> markSet = new HashSet<E>();
		Stack<E> stack = new Stack<E>();
		for (E positiveEdge : graph.getPositiveEdges()){
			if (isBoundaryEdge(positiveEdge)) continue;
			markSet.add(positiveEdge);
			stack.push(positiveEdge);
		}
		while (!stack.isEmpty()) {
			E ab = stack.pop();
			markSet.remove(ab);
			if (!isDelaunay(ab, a)){
				flip(ab, a);
				for (E xy : getPositiveKiteBorder(ab)){
					if (!markSet.contains(xy)){
						markSet.add(xy);
						stack.push(xy);
					}
				}
			}
		}
	}	
	
}
