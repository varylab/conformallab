package de.varylab.discreteconformal.unwrapper.circlepattern;

import static java.lang.Math.exp;

import java.lang.annotation.Annotation;
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
import de.varylab.discreteconformal.unwrapper.quasiisothermic.QuasiisothermicUtility;


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
		HDS extends HalfEdgeDataStructure<V, E, F>,
		GET extends Annotation,
		SET extends Annotation
	> void execute(HDS hds, Map<F, Double> rhoMap, Map<E, Double> thetaMap, AdapterSet a, Class<GET> getter, Class<SET> setter){
		CPEuclideanRotation<V, E, F> rot = new CPEuclideanRotation<V, E, F>();
		calculateGeneric(hds, rot, rhoMap, thetaMap, a, getter, setter);
		// set unlayoutable faces
		List<E> ears = QuasiisothermicUtility.findEarsEdge(hds);
		for (E e : ears) {
			double[] xyFace = a.getD(getter, e.getRightFace());
			a.set(setter, e.getTargetVertex(), xyFace);
		}
	}
	
	static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>,
		GET extends Annotation,
		SET extends Annotation
	> void calculateGeneric(
		HDS hds, 
		Rotation<V, E, F> rot, 
		Map<F, Double> rhoMap, 
		Map<E, Double> thetaMap, 
		AdapterSet a,
		Class<GET> getter, 
		Class<SET> setter
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
		a.set(setter, rootFace, new double[] {0, 0});
		a.set(setter, rootEdge.getTargetVertex(), new double[] {firstPlanarRadius, 0.0});
		layoutEdgeCounterClockwise(firstEdge, rot, rhoMap, thetaMap, a, getter, setter);
		
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
				layoutFace(face, edge, rot, edgeStack, doneEdges, doneFaces, rhoMap, thetaMap, a, getter, setter);
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
		F extends Face<V, E, F>,
		GET extends Annotation,
		SET extends Annotation
	> void layoutFace(
		F face, 
		E edge, 
		Rotation<V, E, F> rot, 
		Stack<E> edgeStack, 
		HashSet<E> doneEdges, 
		HashSet<F> doneFaces,
		Map<F, Double> rhoMap, 
		Map<E, Double> thetaMap, 
		AdapterSet a,
		Class<GET> getter, 
		Class<SET> setter
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
				layoutEdgeCounterClockwise(actEdge, rot, rhoMap, thetaMap, a, getter, setter);
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
				layoutEdgeClockwise(actEdge, rot, rhoMap, thetaMap, a, getter, setter);
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
		F extends Face<V, E, F>,
		GET extends Annotation,
		SET extends Annotation
	> void layoutEdgeClockwise(
		E edge, 
		Rotation<V, E, F> rot, 
		Map<F, Double> rhoMap, 
		Map<E, Double> thetaMap, 
		AdapterSet a,
		Class<GET> getter, 
		Class<SET> setter
	){
		F leftFace = edge.getLeftFace();
		F rightFace = edge.getRightFace();
		V t = edge.getTargetVertex();
		V s = edge.getStartVertex();
		double phi = -rot.getPhi(edge, rhoMap, thetaMap);
		double[] p = a.getD(getter, t);
		double[] c = a.getD(getter, leftFace);
		double[] xy = rot.rotate(p, c, 2*phi, 0.0);
		a.set(setter, s, xy);
		if (rightFace != null){
			double logScale = rhoMap.get(rightFace) - rhoMap.get(leftFace);
			xy = rot.rotate(c, xy, -thetaMap.get(edge), logScale);
			a.set(setter, rightFace, xy);
		}
		
	}
	
	
	/*
	 * Layout startVertex of edge and its right face if non null
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		GET extends Annotation,
		SET extends Annotation
	> void layoutEdgeCounterClockwise(
		E edge, 
		Rotation<V, E, F> rot, 
		Map<F, Double> rhoMap, 
		Map<E, Double> thetaMap, 
		AdapterSet a,
		Class<GET> getter, 
		Class<SET> setter
	){
		F leftFace = edge.getLeftFace();
		F rightFace = edge.getRightFace();
		V t = edge.getTargetVertex();
		V s = edge.getStartVertex();
		double phi = rot.getPhi(edge, rhoMap, thetaMap);
		double[] p = a.getD(getter, s);
		double[] c = a.getD(getter, leftFace);
		double[] xy = rot.rotate(p, c, 2*phi, 0.0);
		a.set(setter, t, xy);
		if (rightFace != null){
			double logScale = rhoMap.get(rightFace) - rhoMap.get(leftFace);
			xy = rot.rotate(c, xy, thetaMap.get(edge), logScale);
			a.set(setter, rightFace, xy);
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
