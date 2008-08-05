package de.varylab.discreteconformal.heds.util;

import java.util.LinkedList;
import java.util.List;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CEdge;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.CVertex;

public class SparseUtility {

	public static int[][] makeNonZeros(CHDS hds) {
		int n = hds.getDomainDimension();
		int[][] nz = new int[n][];
		for (CVertex v : hds.getVertices()) {
			if (v.getSolverIndex() < 0)
				continue;
			List<CEdge> star = HalfEdgeUtils.incomingEdges(v);
			List<Integer> nonZeroIndices = new LinkedList<Integer>();
			nonZeroIndices.add(v.getSolverIndex());
			for (CEdge e : star) {
				CVertex connectedVertex = e.getOppositeEdge().getTargetVertex();
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
