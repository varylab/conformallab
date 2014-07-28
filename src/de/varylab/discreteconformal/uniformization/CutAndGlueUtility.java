package de.varylab.discreteconformal.uniformization;

import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryEdge;

import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;

import de.jreality.math.Matrix;
import de.jreality.math.P2;
import de.jreality.math.Rn;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.CustomVertexInfo;
import de.varylab.discreteconformal.math.RnBig;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class CutAndGlueUtility {

	
	protected static class ValenceComparator implements Comparator<CoVertex> {

		@Override
		public int compare(CoVertex o1, CoVertex o2) {
			int v1 = HalfEdgeUtils.incomingEdges(o1).size();
			int v2 = HalfEdgeUtils.incomingEdges(o2).size();
			return v1 - v2;
		}
		
	}
	
	
	public static Set<CoFace> reglueOutsideFaces(
		CoHDS hds, 
		int maxFaces, 
		FundamentalPolygon p, 
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo,
		int signature
	) {
		int counter = 0;
		Set<CoFace> reglued = new HashSet<CoFace>(); 
		Collection<CoVertex> bList = HalfEdgeUtils.boundaryVertices(hds);
		PriorityQueue<CoVertex> lowValenceFirst = new PriorityQueue<CoVertex>(bList.size(), new ValenceComparator());
		lowValenceFirst.addAll(bList);
		for (CoVertex v : lowValenceFirst) {
			if (!v.isValid()) continue;
			List<CoFace> faces = HalfEdgeUtils.facesIncidentWithVertex(v);
			for (CoFace f : faces) {
				if (counter >= maxFaces) {
					System.out.println(counter + " faces reglued");
					return reglued;
				}
				if (reglued.contains(f)) continue;
				if (HalfEdgeUtils.isInteriorFace(f)) continue;
				if (!isFaceMovable(f)) continue;
				if (!isOutsideFundamentalPolygon(f, p, 0)) continue;
				if (!isFaceMovedToFundamentalDomainByReglue(f, cutInfo, p, 0)) continue;
				try {
					reglueFace(f, cutInfo, signature);
				} catch (AssertionError e) {
					System.err.println(e.getLocalizedMessage());
				}
				reglued.add(f);
				System.out.println("face " + f + " reglued.");
				counter++;
	//			assert HalfEdgeUtils.isValidSurface(hds, true) : "surface should be valid after face reglue";
			}
		}
		System.out.println(counter + " faces reglued");
		return reglued;
	}
	
	public static void reglueSingleFace(CoFace f, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo, int signature) {
		reglueFace(f, cutInfo, signature);
	}
	
	protected static boolean isFaceMovable(CoFace f) {
		for (CoEdge e : HalfEdgeUtils.boundaryEdges(f)) {
			if (HalfEdgeUtils.isBoundaryEdge(e)) continue;
			CoVertex v = e.getNextEdge().getTargetVertex();
			boolean ear = HalfEdgeUtils.incomingEdges(v).size() == 2;
			boolean sb = HalfEdgeUtils.isBoundaryVertex(e.getStartVertex());
			boolean tb = HalfEdgeUtils.isBoundaryVertex(e.getTargetVertex());
			if (sb && tb && !ear) return false;
		}
		return true;
	}
	
	
	
	protected static void reglueFace(
		CoFace f,
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo,
		int signature
	) {
		if (HalfEdgeUtils.isInteriorFace(f)) {
			throw new IllegalArgumentException("can only reglue boundary faces");
		}
		HalfEdgeDataStructure<CoVertex, CoEdge, CoFace> hds = f.getHalfEdgeDataStructure();
		CoEdge e1 = f.getBoundaryEdge();
		CoEdge e2 = e1.getNextEdge();
		CoEdge e3 = e2.getNextEdge();
		int numBoundaryEdges = 0;
		numBoundaryEdges += isBoundaryEdge(e1) ? 1 : 0;
		numBoundaryEdges += isBoundaryEdge(e2) ? 1 : 0;
		numBoundaryEdges += isBoundaryEdge(e3) ? 1 : 0;
		assert numBoundaryEdges < 3 : "a face must have at most two boundary edges";
		if (numBoundaryEdges == 1) {
			CoEdge e = null;
			for (CoEdge be : HalfEdgeUtils.boundaryEdges(f)) {
				if (be.getRightFace() == null) {
					e = be; break;
				}
			}
			assert e != null : "we should have found the boundary edge of face f";
			
			// treat source location
			CoEdge opp = e.getOppositeEdge();
			CoEdge oppNext = opp.getNextEdge();
			CoEdge oppPrev = opp.getPreviousEdge();
			CoEdge eNext = e.getNextEdge();
			CoEdge ePrev = e.getPreviousEdge();
			CoVertex sourceVertex = eNext.getTargetVertex();
			// vertex geometry
			double[] sourcePosT = sourceVertex.T;
			double[] sourcePosP = sourceVertex.P;
			Matrix A = new Matrix(createIsometryFromEdge(e, cutInfo, signature));
			double[] newPos = A.multiplyVector(sourcePosT);
			CoVertex newVertex = hds.addNewVertex();
			newVertex.T = newPos;
			newVertex.P = sourcePosP.clone();
			if (sourceVertex.info != null) {
				newVertex.info = new CustomVertexInfo(sourceVertex.info);
			}
			// relink boundary
			oppPrev.linkNextEdge(eNext);
			ePrev.linkNextEdge(oppNext);
			eNext.setLeftFace(null);
			ePrev.setLeftFace(null);
			// remove old boundary
			hds.removeEdge(e);
			hds.removeEdge(opp);
			// fix incoming edges
			oppPrev.setTargetVertex(eNext.getStartVertex());
			ePrev.setTargetVertex(oppNext.getStartVertex());
			
			// target location
			CoEdge pe = cutInfo.edgeCutMap.get(e);
			CoEdge peOpp = pe.getOppositeEdge();
			CoEdge peOppNext = peOpp.getNextEdge();
			CoEdge peOppPrev = peOpp.getPreviousEdge();
			CoEdge newEdge1 = hds.addNewEdge();
			CoEdge newEdge1Opp = hds.addNewEdge();
			CoEdge newEdge2 = hds.addNewEdge();			
			CoEdge newEdge2Opp = hds.addNewEdge();			
			// link boundary
			peOppPrev.linkNextEdge(newEdge1Opp);
			newEdge1Opp.linkNextEdge(newEdge2Opp);
			newEdge2Opp.linkNextEdge(peOppNext);
			peOpp.linkNextEdge(newEdge2);
			newEdge2.linkNextEdge(newEdge1);
			newEdge1.linkNextEdge(peOpp);
			newEdge2.setTargetVertex(newVertex);
			newEdge1Opp.setTargetVertex(newVertex);
			newEdge1.setTargetVertex(pe.getTargetVertex());
			newEdge2Opp.setTargetVertex(pe.getStartVertex());
			newEdge1.setLeftFace(f);
			newEdge2.setLeftFace(f);
			peOpp.setLeftFace(f);
			newEdge1.linkOppositeEdge(newEdge1Opp);
			newEdge2.linkOppositeEdge(newEdge2Opp);
			
			// fix cut info
			cutInfo.edgeCutMap.put(ePrev, newEdge1Opp);
			cutInfo.edgeCutMap.put(newEdge1Opp, ePrev);
			cutInfo.edgeCutMap.put(eNext, newEdge2Opp);
			cutInfo.edgeCutMap.put(newEdge2Opp, eNext);
			cutInfo.edgeCutMap.put(ePrev.getOppositeEdge(), newEdge1);
			cutInfo.edgeCutMap.put(newEdge1, ePrev.getOppositeEdge());
			cutInfo.edgeCutMap.put(eNext.getOppositeEdge(), newEdge2);
			cutInfo.edgeCutMap.put(newEdge2, eNext.getOppositeEdge());
			cutInfo.edgeCutMap.remove(e);
			cutInfo.edgeCutMap.remove(opp);
			cutInfo.edgeCutMap.remove(pe);
			cutInfo.edgeCutMap.remove(peOpp);
			cutInfo.vertexCopyMap.put(sourceVertex, newVertex);
		} else {
			CoEdge e = null;
			for (CoEdge be : HalfEdgeUtils.boundaryEdges(f)) {
				if (be.getRightFace() != null) {
					e = be; break;
				}
			}
			assert e != null : "we should have found the non-boundary edge of face f";
			CoEdge b1 = e.getNextEdge();
			CoEdge b2 = e.getPreviousEdge();
			CoEdge pb1 = cutInfo.edgeCutMap.get(b1);
			CoEdge pb2 = cutInfo.edgeCutMap.get(b2);
			if (pb1.getOppositeEdge().getNextEdge() != pb2.getOppositeEdge()) {
				System.out.println("no continuous edge identification at face " + f);
				return;
			}
			CoEdge b1Opp = b1.getOppositeEdge();
			CoEdge b2Opp = b2.getOppositeEdge();
			CoEdge b1OppNext = b1Opp.getNextEdge();
			CoEdge b2OppPrev = b2Opp.getPreviousEdge();
			CoVertex v = b1.getTargetVertex();
			CoVertex v1 = b1.getStartVertex();
			CoVertex v2 = e.getStartVertex();
			// relink source boundary
			b2OppPrev.linkNextEdge(e);
			e.linkNextEdge(b1OppNext);
			e.setLeftFace(null);
			// delete old boundary
			hds.removeEdge(b1);
			hds.removeEdge(b2);
			hds.removeEdge(b1Opp);
			hds.removeEdge(b2Opp);
			hds.removeVertex(v);
			// fix incoming edges
			b2OppPrev.setTargetVertex(v2);
			e.setTargetVertex(v1);
			
			// target location
			CoEdge newEdge = hds.addNewEdge();
			CoEdge newEdgeOpp = hds.addNewEdge();
			CoEdge pb1Opp = pb1.getOppositeEdge();
			CoEdge pb2Opp = pb2.getOppositeEdge();
			CoEdge pb1OppPrev = pb1Opp.getPreviousEdge();
			CoEdge pb2OppNext = pb2Opp.getNextEdge();
			// link new boundary
			pb1OppPrev.linkNextEdge(newEdgeOpp);
			newEdgeOpp.linkNextEdge(pb2OppNext);
			// relink face
			newEdge.linkOppositeEdge(newEdgeOpp);
			pb2Opp.linkNextEdge(newEdge);
			newEdge.linkNextEdge(pb1Opp);
			newEdge.setLeftFace(f);
			pb1Opp.setLeftFace(f);
			pb2Opp.setLeftFace(f);
			newEdge.setTargetVertex(pb1.getTargetVertex());
			newEdgeOpp.setTargetVertex(pb2.getStartVertex());
			
			// fix cut info
			cutInfo.edgeCutMap.put(e, newEdgeOpp);
			cutInfo.edgeCutMap.put(newEdgeOpp, e);
			cutInfo.edgeCutMap.put(e.getOppositeEdge(), newEdge);
			cutInfo.edgeCutMap.put(newEdge, e.getOppositeEdge());
			cutInfo.edgeCutMap.remove(b1);
			cutInfo.edgeCutMap.remove(b2);
			cutInfo.edgeCutMap.remove(b1Opp);
			cutInfo.edgeCutMap.remove(b1Opp);
			cutInfo.vertexCopyMap.remove(v);
		}
	}
	

	protected static boolean isOutsideFundamentalPolygon(CoFace f, FundamentalPolygon p, double tol) {
		for (CoVertex v : HalfEdgeUtils.boundaryVertices(f)) {
			if (isInsideFundamentalPolygon(v, p, tol)) return false;
		}
		return true;
	}
	
	protected static boolean isInsideFundamentalPolygon(CoVertex v, FundamentalPolygon p, double tol) {
		double[] vt = P2.projectP3ToP2(null, v.T); 
		for (FundamentalEdge e : p.getEdges()) {
			double[] s = P2.projectP3ToP2(null, RnBig.toDouble(null, e.startPosition));
			double[] t = P2.projectP3ToP2(null, RnBig.toDouble(null, e.nextEdge.startPosition));
			double[] line = P2.lineFromPoints(null, s, t);
			double dot = Rn.innerProduct(line, vt);
			if (dot < tol) return false;
		}
		return true;
	}
	
	protected static boolean isFaceMovedToFundamentalDomainByReglue(
		CoFace f, 
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo, 
		FundamentalPolygon p, 
		double tol
	) {
		for (CoEdge e : HalfEdgeUtils.boundaryEdges(f)) {
			CoEdge pe = cutInfo.edgeCutMap.get(e);
			if (pe != null) {
				if (!isInsideFundamentalPolygon(pe.getStartVertex(), p, tol)) return false;
				if (!isInsideFundamentalPolygon(pe.getTargetVertex(), p, tol)) return false;
			}
		}
		return true;
	}
	
	
	protected static double[] createIsometryFromEdge(CoEdge e, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo, int signature) {
		if (!HalfEdgeUtils.isBoundaryEdge(e)) {
			throw new IllegalArgumentException("No boundary edge");
		}
		CoEdge ePartner = cutInfo.edgeCutMap.get(e);
		assert ePartner != null : "every boundary edge has to be in the cut info";
		double[] s1 = P2.projectP3ToP2(null, e.getStartVertex().T);
		double[] t1 = P2.projectP3ToP2(null, e.getTargetVertex().T);
		double[] s2 = P2.projectP3ToP2(null, ePartner.getStartVertex().T);
		double[] t2 = P2.projectP3ToP2(null, ePartner.getTargetVertex().T);
		double[] a = P2.makeDirectIsometryFromFrames(null, s1, t1, t2, s2, signature);
		return P2.imbedMatrixP2InP3(null, a);
	}

	
}
