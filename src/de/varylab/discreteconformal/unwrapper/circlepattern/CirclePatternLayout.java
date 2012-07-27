package de.varylab.discreteconformal.unwrapper.circlepattern;

import static java.lang.Math.exp;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.generic.Position2d;
import de.varylab.discreteconformal.unwrapper.isothermic.IsothermicUtility;


/**
 * Copyright 2005 <a href="http://www.sechel.de">Stefan Sechelmann</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Stefan Sechelmann
 * @see koebe.KoebePolyhedron
 */
public final class CirclePatternLayout {

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds, Map<F, Double> rhoMap, Map<E, Double> thetaMap, AdapterSet a){
		CPEuclideanRotation<V, E, F> rot = new CPEuclideanRotation<V, E, F>();
		calculateGeneric(hds, rot, rhoMap, thetaMap, a);
		// set unlayoutable faces
		List<E> ears = IsothermicUtility.findEarsEdge(hds);
		for (E e : ears) {
			double[] xyFace = a.getD(Position2d.class, e.getRightFace());
			a.set(Position.class, e.getTargetVertex(), xyFace);
		}
	}
	
	static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void calculateGeneric(
		HDS hds, 
		Rotation<V, E, F> rot, 
		Map<F, Double> rhoMap, 
		Map<E, Double> thetaMap, 
		AdapterSet a
	) { 
		Stack<E> edgeStack = new Stack<E>();
		HashSet<E> doneEdges = new HashSet<E>();
		HashSet<F> doneFaces = new HashSet<F>();	
		
		// Init ---------------------------------------
		F rootFace = hds.getFace(0);
		for (F f : hds.getFaces()) {
			if (HalfEdgeUtils.isInteriorFace(f)){
				rootFace = f;
				break;
			}
		}
		E rootEdge = rootFace.getBoundaryEdge();
		E firstEdge = rootEdge.getNextEdge();

		double firstPlanarRadius = exp(rhoMap.get(rootFace));
		a.set(Position.class, rootFace, new double[] {0, 0});
		a.set(Position.class, rootEdge.getTargetVertex(), new double[] {firstPlanarRadius, 0.0});
		layoutEdgeCounterClockwise(firstEdge, rot, rhoMap, thetaMap, a);
		
		if (firstEdge.getRightFace() != null)
			edgeStack.push(firstEdge.getOppositeEdge());
		edgeStack.push(firstEdge);
		
		doneEdges.add(firstEdge);
		doneEdges.add(firstEdge.getOppositeEdge());
		// ---------------------------------------------
		
		
		while (!edgeStack.isEmpty()){
			E edge = edgeStack.pop();
			F face = edge.getLeftFace();
			if (!doneFaces.contains(face)) {
				layoutFace(face, edge, rot, edgeStack, doneEdges, doneFaces, rhoMap, thetaMap, a);
			}
		}
	}
	
	/*
	 * Layouts all edges and surrounding faces of face
	 * edge and face must have been layouted already
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void layoutFace(
		F face, 
		E edge, 
		Rotation<V, E, F> rot, 
		Stack<E> edgeStack, 
		HashSet<E> doneEdges, 
		HashSet<F> doneFaces,
		Map<F, Double> rhoMap, 
		Map<E, Double> thetaMap, 
		AdapterSet a
	){
		doneFaces.add(face);
		boolean stoppedAtBoundary = false;
		E actEdge = edge.getNextEdge();
		// layout clockwise
		while (actEdge != edge) {
			if (!HalfEdgeUtils.isInteriorEdge(actEdge)) {
				stoppedAtBoundary = true;
				break;
			}
			if (!doneEdges.contains(actEdge)){
				layoutEdgeCounterClockwise(actEdge, rot, rhoMap, thetaMap, a);
				if (actEdge.getRightFace() != null)
					edgeStack.push(actEdge.getOppositeEdge());
				doneEdges.add(actEdge);
				doneEdges.add(actEdge.getOppositeEdge());
			}
			actEdge = actEdge.getNextEdge();
		}
		if (!stoppedAtBoundary)
			return;
		// layout counter clockwise if we need to
		actEdge = edge.getPreviousEdge();
		while (actEdge != edge) {
			if (!HalfEdgeUtils.isInteriorEdge(actEdge))
				return;
			if (!doneEdges.contains(actEdge)){
				layoutEdgeClockwise(actEdge, rot, rhoMap, thetaMap, a);
				if (actEdge.getRightFace() != null)
					edgeStack.push(actEdge.getOppositeEdge());
				doneEdges.add(actEdge);
				doneEdges.add(actEdge.getOppositeEdge());
			}
			actEdge = actEdge.getPreviousEdge();
		}
	}
	
	
	
	/*
	 * Layout startVertex of edge and its right face if non null
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void layoutEdgeClockwise(
		E edge, 
		Rotation<V, E, F> rot, 
		Map<F, Double> rhoMap, 
		Map<E, Double> thetaMap, 
		AdapterSet a
	){
		F leftFace = edge.getLeftFace();
		F rightFace = edge.getRightFace();
		V t = edge.getTargetVertex();
		V s = edge.getStartVertex();
		double phi = -rot.getPhi(edge, rhoMap, thetaMap);
		double[] p = a.getD(Position2d.class, t);
		double[] c = a.getD(Position2d.class, leftFace);
		double[] xy = rot.rotate(p, c, 2*phi, 0.0);
		a.set(Position.class, s, xy);
		if (rightFace != null){
			double logScale = rhoMap.get(rightFace) - rhoMap.get(leftFace);
			xy = rot.rotate(c, xy, -thetaMap.get(edge), logScale);
			a.set(Position.class, rightFace, xy);
		}
		
	}
	
	
	/*
	 * Layout startVertex of edge and its right face if non null
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void layoutEdgeCounterClockwise(
		E edge, 
		Rotation<V, E, F> rot, 
		Map<F, Double> rhoMap, 
		Map<E, Double> thetaMap, 
		AdapterSet a
	){
		F leftFace = edge.getLeftFace();
		F rightFace = edge.getRightFace();
		V t = edge.getTargetVertex();
		V s = edge.getStartVertex();
		double phi = rot.getPhi(edge, rhoMap, thetaMap);
		double[] p = a.getD(Position2d.class, s);
		double[] c = a.getD(Position2d.class, leftFace);
		double[] xy = rot.rotate(p, c, 2*phi, 0.0);
		a.set(Position.class, t, xy);
		if (rightFace != null){
			double logScale = rhoMap.get(rightFace) - rhoMap.get(leftFace);
			xy = rot.rotate(c, xy, thetaMap.get(edge), logScale);
			a.set(Position.class, rightFace, xy);
		}
		
	}
	
	protected static interface Rotation <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> {
		public double[] rotate(double[] p, double[] center, Double phi, Double logScale);
		public double getPhi(E edge, Map<F, Double> rhoMap, Map<E, Double> thetaMap);
		public Double getRadius(Double rho);
	}
	
}
