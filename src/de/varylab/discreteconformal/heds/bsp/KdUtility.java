package de.varylab.discreteconformal.heds.bsp;

import static de.jtem.halfedge.util.HalfEdgeUtils.facesIncidentWithVertex;
import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;

import java.util.HashSet;
import java.util.Vector;

import de.varylab.discreteconformal.heds.CEdge;
import de.varylab.discreteconformal.heds.CFace;
import de.varylab.discreteconformal.heds.CVertex;

public final class KdUtility {

	public static Vector<CFace> collectFacesInRadius(KdTree<CVertex> kdTree, HasPosition p, double radius) {
		Vector<CVertex> vertresult = kdTree.collectInRadius(p, radius);
		HashSet<CFace> faceSet = new HashSet<CFace>(vertresult.size() * 2);
		
		for (CVertex v : vertresult)
			faceSet.addAll(facesIncidentWithVertex(v));
		
		Vector<CFace> result = new Vector<CFace>();
		result.addAll(faceSet);
		return result;
	}
	
	
	public static Vector<CEdge> collectEdgesInRadius(KdTree<CVertex> kdTree, HasPosition p, double radius) {
		Vector<CVertex> vertresult = kdTree.collectInRadius(p, radius);
		HashSet<CEdge> edgeSet = new HashSet<CEdge>(vertresult.size() * 3);
		
		for (CVertex v : vertresult)
			edgeSet.addAll(incomingEdges(v));
		
		Vector<CEdge> result = new Vector<CEdge>();
		result.addAll(edgeSet);
		return result;
	}
	
}
