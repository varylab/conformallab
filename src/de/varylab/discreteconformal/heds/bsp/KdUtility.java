package de.varylab.discreteconformal.heds.bsp;

import static de.jtem.halfedge.util.HalfEdgeUtils.facesIncidentWithVertex;
import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;

import java.util.HashSet;
import java.util.Vector;

import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

public final class KdUtility {

	public static Vector<CoFace> collectFacesInRadius(KdTree<CoVertex> kdTree, HasPosition p, double radius) {
		Vector<CoVertex> vertresult = kdTree.collectInRadius(p, radius);
		HashSet<CoFace> faceSet = new HashSet<CoFace>(vertresult.size() * 2);
		
		for (CoVertex v : vertresult)
			faceSet.addAll(facesIncidentWithVertex(v));
		
		Vector<CoFace> result = new Vector<CoFace>();
		result.addAll(faceSet);
		return result;
	}
	
	
	public static Vector<CoEdge> collectEdgesInRadius(KdTree<CoVertex> kdTree, HasPosition p, double radius) {
		Vector<CoVertex> vertresult = kdTree.collectInRadius(p, radius);
		HashSet<CoEdge> edgeSet = new HashSet<CoEdge>(vertresult.size() * 3);
		
		for (CoVertex v : vertresult)
			edgeSet.addAll(incomingEdges(v));
		
		Vector<CoEdge> result = new Vector<CoEdge>();
		result.addAll(edgeSet);
		return result;
	}
	
}
