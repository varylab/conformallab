package de.varylab.discreteconformal.heds;

import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryVertex;
import static java.lang.Math.abs;

import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.CG;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import de.varylab.discreteconformal.heds.util.GraphUtility;
import de.varylab.discreteconformal.heds.util.Search;

public class CCones {

	
	/**
	 * Cuts the mesh along paths from the vertices of cones to the boundary
	 * @param hds the mesh
	 * @param cones the cones
	 */
	public static void cutMesh(CHDS hds, Collection<CVertex> cones) {
		Set<CVertex> bSet = new HashSet<CVertex>();
		for (CVertex v : hds.getVertices()){
			if (isBoundaryVertex(v))
				bSet.add(v);
		}

		LambdaEdgeComparatore comp = new LambdaEdgeComparatore();
		for (CVertex c : cones) {
			//TODO Implement shortest path cutting
			java.util.Vector<CEdge> path = Search.bFS(c, bSet, true, comp);
			for (CEdge e : path) {
				CEdge eOpp = e.getOppositeEdge();
				Map<CVertex, CVertex> vMap = GraphUtility.cutAtEdge(e);
				for (CVertex v : vMap.keySet()) {
					CVertex nV = vMap.get(v);
					nV.setSolverIndex(v.getSolverIndex());
					nV.setTheta(v.getTheta());
					nV.setPosition(v.getPosition());
				}
				e.getOppositeEdge().setLambda(e.getLambda());
				eOpp.getOppositeEdge().setLambda(eOpp.getLambda());
			}
		}
	}
	
	
	protected static class LambdaEdgeComparatore implements Comparator<CEdge> {

		public int compare(CEdge o1, CEdge o2) {
			return o1.getLambda() < o2.getLambda() ? -1 : 1;
		}
		
	}
	
	
	/**
	 * calculates the possibly best cone points in hds and sets up new solver indices.
	 * Needs an invariant data prepared CHDS
	 * @param hds
	 * @param cones
	 */
	public static Collection<CVertex> setUpMesh(CHDS hds, int cones) {
		Collection<CVertex> result = new LinkedList<CVertex>();
		for (int i = 0; i < cones; i++) {
			int n = hds.getDomainDimension();
			Vector u = new DenseVector(n);
			Vector G = new DenseVector(n);
			Matrix H = new DenseMatrix(n, n);
			CG cg = new CG(u);
			hds.conformalEnergy(u, null, G, H);
			try {
				cg.solve(H, G, u);
			} catch (IterativeSolverNotConvergedException e) {
				e.printStackTrace();
				return result;
			}
			int max = getMaxAbsIndex(u);
			CVertex coneVertex = findVertexWidthSolverIndex(hds, max);
			coneVertex.setSolverIndex(-1);
			reorderSolverIndices(hds);
			result.add(coneVertex);
		}
		return result;
	}
	
	
	/**
	 * Reorders the solver indices for  
	 * @param hds
	 */
	public static void reorderSolverIndices(CHDS hds) {
		int i = 0;
		for (CVertex v : hds.getVertices()) {
			if (v.getSolverIndex() >= 0)
				v.setSolverIndex(i++);
		}
		hds.setDomainDimension(hds.getDomainDimension() - 1);
	}
	
	
	
	
	/**
	 * returns the index of the element with the largest absolute value in u
	 * @param u
	 * @return
	 */
	public static int getMaxAbsIndex(Vector u) {
		int max = 0;
		double maxVal = abs(u.get(max));
		for (int i = 1; i < u.size(); i++) {
			double val = abs(u.get(i));
			if (maxVal < val) {
				max = i;
				maxVal = val;
			}
		}
		return max;
	}
	
	
	/**
	 * Returns the first vertex with solver index i
	 * @param hds
	 * @param i
	 * @return
	 */
	public static CVertex findVertexWidthSolverIndex(CHDS hds, int i) {
		for (CVertex v : hds.getVertices()) {
			if (v.getSolverIndex() == i)
				return v;
		}
		return null;
	}
	
	
}
