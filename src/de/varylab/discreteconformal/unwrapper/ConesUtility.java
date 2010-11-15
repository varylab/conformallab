package de.varylab.discreteconformal.unwrapper;

import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryVertex;
import static de.varylab.discreteconformal.unwrapper.EuclideanLayout.calculateAngleSum;
import static de.varylab.discreteconformal.util.SparseUtility.makeNonZeros;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.exp;

import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.CG;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanOptimizable;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.Search;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

public class ConesUtility {

	
	/**
	 * Cuts the mesh along paths from the vertices of cones to the boundary
	 * @param hds the mesh
	 * @param cones the cones
	 */
	public static void cutMesh(CoHDS hds, Collection<CoVertex> cones, Vector u) {
		if (cones.isEmpty()) return;

		for (CoVertex c : cones) {
			// get boundary set
			Set<CoVertex> bSet = new HashSet<CoVertex>();
			for (CoVertex v : hds.getVertices()){
				if (isBoundaryVertex(v)) {
					bSet.add(v);
				}
			}
			// find path to boundary
			List<CoEdge> path = Search.getShortestPath(c, bSet, new EdgeLengthAdapter(u));
			// cut
			for (CoEdge e : path) {
				CoEdge eOpp = e.getOppositeEdge();
				Map<CoVertex, CoVertex> vMap = CuttingUtility.cutAtEdge(e);
				for (CoVertex v : vMap.keySet()) {
					CoVertex nV = vMap.get(v);
					nV.setSolverIndex(v.getSolverIndex());
					nV.setTheta(v.getTheta());
					nV.P = v.P.clone();
				}
				e.getOppositeEdge().setLambda(e.getLambda());
				eOpp.getOppositeEdge().setLambda(eOpp.getLambda());
			}
		}

	}
	
	
	protected static class EdgeLengthAdapter implements WeightAdapter<CoEdge> {

		private Vector
			u = null;
		
		public EdgeLengthAdapter(Vector u) {
			this.u = u;
		}
		
		@Override
		public double getWeight(CoEdge e) {
			CoVertex v1 = e.getStartVertex();
			CoVertex v2 = e.getTargetVertex();
			Double u1 = v1.getSolverIndex() >= 0 ? u.get(v1.getSolverIndex()) : 0.0; 
			Double u2 = v2.getSolverIndex() >= 0 ? u.get(v2.getSolverIndex()) : 0.0;
			Double lambda = e.getLambda();
			return exp(lambda + u1 + u2);
		}
		
	}
	
	
	protected static class LambdaEdgeComparatore implements Comparator<CoEdge> {

		@Override
		public int compare(CoEdge o1, CoEdge o2) {
			return o1.getLambda() < o2.getLambda() ? -1 : 1;
		}
		
	}
	
	
	/**
	 * calculates the possibly best cone points in hds and sets up new solver indices.
	 * Needs an invariant data prepared CHDS
	 * @param hds
	 * @param numAutomaticCones
	 */
	public static Collection<CoVertex> setUpCones(CoHDS hds, int numAutomaticCones) {
		Collection<CoVertex> result = new LinkedList<CoVertex>();
		// automatic cones
		for (int i = 0; i < numAutomaticCones; i++) {
			CEuclideanOptimizable opt = new CEuclideanOptimizable(hds);
			int n = opt.getDomainDimension();
			Vector u = new DenseVector(n);
			Vector G = new DenseVector(n);
			Matrix H = new CompRowMatrix(n, n, makeNonZeros(hds));
			CG cg = new CG(u);
			opt.evaluate(u, G, H);
			try {
				cg.solve(H, G, u);
			} catch (IterativeSolverNotConvergedException e) {
				e.printStackTrace();
				return result;
			}
			int max = getMaxAbsIndex(u);
			CoVertex coneVertex = findVertexWidthSolverIndex(hds, max);
			coneVertex.setSolverIndex(-1);
			coneVertex.setTheta(0.0);
			result.add(coneVertex);
			reorderSolverIndices(hds);
		}
		// custom cones
		for (CoVertex v : hds.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) continue;
			if (v.info == null) continue;
			if (v.info.useCustomTheta) {
				v.setSolverIndex(0);
				v.setTheta(v.info.theta);
				result.add(v);
			}
		}
		reorderSolverIndices(hds);
		return result;
	}
	
	
	/**
	 * Reorders the solver indices for  
	 * @param hds
	 */
	public static void reorderSolverIndices(CoHDS hds) {
		int i = 0;
		for (CoVertex v : hds.getVertices()) {
			if (v.getSolverIndex() >= 0) {
				v.setSolverIndex(i++);
			}
		}
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
	public static CoVertex findVertexWidthSolverIndex(CoHDS hds, int i) {
		for (CoVertex v : hds.getVertices()) {
			if (v.getSolverIndex() == i)
				return v;
		}
		return null;
	}
	
	
	public static Collection<CoVertex> quantizeCones(CoHDS hds, Collection<CoVertex> cones, QuantizationMode qm) {
		List<CoVertex> result = new LinkedList<CoVertex>(cones);
		for (CoVertex v : cones) {
			double theta = calculateAngleSum(v);
			switch (qm) {
			case AllAngles:
				break;
			case Straight:
				theta = PI;
				break;
			case QuadsStrict:
				if (theta < PI / 4) {
					theta = 2 * PI;
					result.remove(v);
				} else {
					double q = theta % (PI/2);
					if (q > PI/4) {
						theta += (PI/2 - q);
					} else {
						theta -= q;
					}
				}
			case Quads:	
				if (theta < PI / 8) {
					theta = 2 * PI;
					result.remove(v);
				} else {
					double q = theta % (PI/4);
					if (q > PI/8) {
						theta += (PI/4 - q);
					} else {
						theta -= q;
					}
				}
				break;
			case Hexagons:
				if (theta < PI / 12) {
					theta = 2 * PI;
					result.remove(v);
				} else {
					double q = theta % (PI/6);
					if (q > PI/12) {
						theta += (PI/6 - q);
					} else {
						theta -= q;
					}
				}
				break;
			}
			v.setTheta(theta);
			v.setSolverIndex(0);
		}
		reorderSolverIndices(hds);
		return result;
	}
	
	
}
