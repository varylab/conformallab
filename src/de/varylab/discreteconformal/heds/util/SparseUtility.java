package de.varylab.discreteconformal.heds.util;

import java.util.LinkedList;
import java.util.List;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class SparseUtility {

	
	
	public static int[] getPETScNonZeros(CoHDS hds){
		int [][] sparseStucture = makeNonZeros(hds);
		int [] nnz = new int[sparseStucture.length];
		for(int i = 0; i < nnz.length; i++){
			nnz[i] = sparseStucture[i].length;
		}
		return nnz;
	}
	
	
	public static int[][] makeNonZeros(CoHDS hds) {
		int n = 0;
		for (CoVertex v : hds.getVertices()) {
			if (v.getSolverIndex() >= 0) {
				n++;
			}
		}
		int[][] nz = new int[n][];
		for (CoVertex v : hds.getVertices()) {
			if (v.getSolverIndex() < 0)
				continue;
			List<CoEdge> star = HalfEdgeUtils.incomingEdges(v);
			List<Integer> nonZeroIndices = new LinkedList<Integer>();
			nonZeroIndices.add(v.getSolverIndex());
			for (CoEdge e : star) {
				CoVertex connectedVertex = e.getOppositeEdge().getTargetVertex();
				if (connectedVertex.getSolverIndex() < 0)
					continue;
				nonZeroIndices.add(connectedVertex.getSolverIndex());
			}
			nz[v.getSolverIndex()] = new int[nonZeroIndices.size()];
			int i = 0;
			for (Integer index : nonZeroIndices) {
				nz[v.getSolverIndex()][i++] = index;
			}
		}
		return nz;
	}

}
