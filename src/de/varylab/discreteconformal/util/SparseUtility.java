package de.varylab.discreteconformal.util;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class SparseUtility {

	
	public static int[] getPETScNonZeros(int[][] sparseStucture){
		int [] nnz = new int[sparseStucture.length];
		for(int i = 0; i < nnz.length; i++){
			nnz[i] = sparseStucture[i].length;
		}
		return nnz;
	}
	
	
	public static int[] getPETScNonZeros(CoHDS hds){
		int[][] sparseStucture = makeNonZeros(hds);
		int[] nnz = new int[sparseStucture.length];
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
		for (CoEdge e : hds.getPositiveEdges()) {
			if (e.getSolverIndex() >= 0) {
				n++;
			}
		}
		Map<Integer, TreeSet<Integer>> nonZeros = new HashMap<Integer, TreeSet<Integer>>();
		for (int i = 0; i < n; i++) {
			nonZeros.put(i, new TreeSet<Integer>());
		}
		for (CoVertex v : hds.getVertices()) {
			if (v.getSolverIndex() < 0) {
				continue;
			}
			List<CoEdge> star = HalfEdgeUtils.incomingEdges(v);
			Set<Integer> nonZeroIndices = nonZeros.get(v.getSolverIndex());
			nonZeroIndices.add(v.getSolverIndex());
			for (CoEdge e : star) {
				CoVertex connectedVertex = e.getOppositeEdge().getTargetVertex();
				if (connectedVertex.getSolverIndex() >= 0) {
					nonZeroIndices.add(connectedVertex.getSolverIndex());
				}
				if (e.getSolverIndex() >= 0) {
					nonZeroIndices.add(e.getSolverIndex());
				}
				if (e.getPreviousEdge().getSolverIndex() >= 0) {
					nonZeroIndices.add(e.getPreviousEdge().getSolverIndex());
				}
			}
		}
		for (CoEdge e : hds.getEdges()) {
			if (e.getSolverIndex() < 0) {
				continue;
			}
			Set<Integer> nonZeroIndices = nonZeros.get(e.getSolverIndex());
			
			// quadratic derivative
			nonZeroIndices.add(e.getSolverIndex());
			
			// mixed edge derivatives
			if (e.getNextEdge().getSolverIndex() >= 0) {
				nonZeroIndices.add(e.getNextEdge().getSolverIndex());
			}
			if (e.getPreviousEdge().getSolverIndex() >= 0) {
				nonZeroIndices.add(e.getPreviousEdge().getSolverIndex());
			}
			
			// mixed vertex derivatives
			if (e.getTargetVertex().getSolverIndex() >= 0) {
				nonZeroIndices.add(e.getTargetVertex().getSolverIndex());
			}
			if (e.getNextEdge().getTargetVertex().getSolverIndex() >= 0) {
				nonZeroIndices.add(e.getNextEdge().getTargetVertex().getSolverIndex());
			}
		}
		int[][] nz = new int[n][];
		for (int j = 0; j < n; j++) {
			Set<Integer> nonZeroIndices = nonZeros.get(j);
			nz[j] = new int[nonZeroIndices.size()];
			int i = 0;
			for (Integer index : nonZeroIndices) {
				nz[j][i++] = index;
			}
		}
		return nz;
	}

}
