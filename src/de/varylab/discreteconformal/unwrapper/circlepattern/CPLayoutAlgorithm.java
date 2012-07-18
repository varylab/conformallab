/**
This file is part of a jTEM project.
All jTEM projects are licensed under the FreeBSD license 
or 2-clause BSD license (see http://www.opensource.org/licenses/bsd-license.php). 

Copyright (c) 2002-2010, Technische Universität Berlin, jTEM
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

-	Redistributions of source code must retain the above copyright notice, 
	this list of conditions and the following disclaimer.

-	Redistributions in binary form must reproduce the above copyright notice, 
	this list of conditions and the following disclaimer in the documentation 
	and/or other materials provided with the distribution.
 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS 
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY 
OF SUCH DAMAGE.
**/

package de.varylab.discreteconformal.unwrapper.circlepattern;

import static java.lang.Math.exp;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import javax.vecmath.Point2d;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.unwrapper.circlepattern.CPLayoutAdapters.XYFace;
import de.varylab.discreteconformal.unwrapper.circlepattern.CPLayoutAdapters.XYVertex;
import de.varylab.discreteconformal.unwrapper.isothermic.IsothermicUtility;


/**
 * A layouter for circle patterns calculated with 
 * koebe.KoebePolyhedron
 * <p>
 * Copyright 2005 <a href="http://www.sechel.de">Stefan Sechelmann</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Stefan Sechelmann
 * @see koebe.KoebePolyhedron
 */
public class CPLayoutAlgorithm <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>
> {

	private XYVertex<V> 
		xyVertex = null;
	private XYFace<F> 
		xyFace = null;
	private Map<F, Double> 
		rhoMap = null;
	private Map<E, Double> 
		thetaMap = null;
	
	
	public CPLayoutAlgorithm(XYVertex<V> xyV, XYFace<F> xyF, Map<F, Double> rhoMap, Map<E, Double> thetaMap) {
		this.xyVertex = xyV;
		this.xyFace = xyF;
		this.thetaMap = thetaMap;
		this.rhoMap = rhoMap;
	}
	

	public <
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(
		 HDS hds
	){
		CPEuclideanRotation<V, E, F> rot = new CPEuclideanRotation<V, E, F>();
		calculateGeneric(hds, rot);
		// set unlayoutable faces
		List<E> ears = IsothermicUtility.findEarsEdge(hds);
		for (E e : ears) {
			Point2d xy = xyFace.getXY(e.getRightFace(), new Point2d());
			xyVertex.setXY(e.getTargetVertex(), xy);
		}
	}
	
	
	private <
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void calculateGeneric(
			HDS hds, 
			Rotation<V, E, F> rot
	){
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
		
		xyFace.setXY(rootFace, new Point2d());
		double firstPlanarRadius = exp(rhoMap.get(rootFace));
		xyVertex.setXY(rootEdge.getTargetVertex(), new Point2d(firstPlanarRadius, 0.0));
		layoutEdgeCounterClockwise(firstEdge, rot);
		
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
				layoutFace(face, edge, rot, edgeStack, doneEdges, doneFaces);
			}
		}
	}
	
	/*
	 * Layouts all edges and surrounding faces of face
	 * edge and face must have been layouted already
	 */
	private void layoutFace(
		F face, 
		E edge, 
		Rotation<V, E, F> rot, 
		Stack<E> edgeStack, 
		HashSet<E> doneEdges, 
		HashSet<F> doneFaces
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
				layoutEdgeCounterClockwise(actEdge, rot);
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
				layoutEdgeClockwise(actEdge, rot);
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
	private void layoutEdgeClockwise(
		E edge, 
		Rotation<V, E, F> rot
	){
		F leftFace = edge.getLeftFace();
		F rightFace = edge.getRightFace();
		V t = edge.getTargetVertex();
		V s = edge.getStartVertex();
		Double phi = -rot.getPhi(edge, rhoMap, thetaMap);
		Point2d xy = rot.rotate(xyVertex.getXY(t, new Point2d()), xyFace.getXY(leftFace, new Point2d()), 2*phi, 0.0);
		xyVertex.setXY(s, xy);	
		if (rightFace != null){
			Double logScale = rhoMap.get(rightFace) - rhoMap.get(leftFace);
			xy = rot.rotate(xyFace.getXY(leftFace, new Point2d()), xyVertex.getXY(s, new Point2d()), -thetaMap.get(edge), logScale);
			xyFace.setXY(rightFace, xy);
		}
		
	}
	
	
	/*
	 * Layout startVertex of edge and its right face if non null
	 */
	private void layoutEdgeCounterClockwise(
		E edge, 
		Rotation<V, E, F> rot
	){
		F leftFace = edge.getLeftFace();
		F rightFace = edge.getRightFace();
		V t = edge.getTargetVertex();
		V s = edge.getStartVertex();
		Double phi = rot.getPhi(edge, rhoMap, thetaMap);
		Point2d xy = rot.rotate(xyVertex.getXY(s, new Point2d()), xyFace.getXY(leftFace, new Point2d()), 2*phi, 0.0);
		xyVertex.setXY(t, xy);	
		if (rightFace != null){
			Double logScale = rhoMap.get(rightFace) - rhoMap.get(leftFace);
			xyFace.setXY(rightFace, rot.rotate(xyFace.getXY(leftFace, new Point2d()), xyVertex.getXY(t, new Point2d()), thetaMap.get(edge), logScale));
		}
		
	}
	
	protected static interface Rotation <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> {
		public Point2d rotate(Point2d p, Point2d center, Double phi, Double logScale);
		public double getPhi(E edge, Map<F, Double> rhoMap, Map<E, Double> thetaMap);
		public Double getRadius(Double rho);
	}
	
}