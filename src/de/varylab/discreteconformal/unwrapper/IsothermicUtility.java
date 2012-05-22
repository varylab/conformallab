package de.varylab.discreteconformal.unwrapper;

import static java.lang.Math.PI;
import static java.lang.Math.abs;

import java.util.HashMap;
import java.util.Map;

import javax.vecmath.Point2d;

import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.functional.Functional;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.PETSc;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.jtem.jtao.Tao.Method;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.CPLayoutAdapters.XYFace;
import de.varylab.discreteconformal.unwrapper.CPLayoutAdapters.XYVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CPEuclideanApplication;

public class IsothermicUtility {

	public static class CPLayoutAdapters implements XYVertex<CoVertex>, XYFace<CoFace> {
	
		@Override
		public Point2d getXY(CoFace v, Point2d xy) {
			return new Point2d(new double[] {0,0});
		}
	
		@Override
		public void setXY(CoFace f, Point2d xy) {
		}
	
		@Override
		public Point2d getXY(CoVertex v, Point2d xy) {
			return new Point2d(new double[] {v.T[0], v.T[1]});
		}
	
		@Override
		public void setXY(CoVertex v, Point2d xy) {
			v.T[0] = xy.x;
			v.T[1] = xy.y;
			v.T[2] = 0.0;
			v.T[3] = 1.0;
		}
		
	}

	/**
	 * Calculate the angle between the edges that belong to alpha1 and alpha2.
	 * @param alpha1
	 * @param alpha2
	 * @param alpha3
	 * @return
	 */
	public static double calculateTriangleAngle(double alpha1, double alpha2, double alpha3) {
		alpha1 = IsothermicUtility.normalizeAngle(alpha1);
		alpha2 = IsothermicUtility.normalizeAngle(alpha2);
		alpha3 = IsothermicUtility.normalizeAngle(alpha3);
		double beta = abs(alpha2 - alpha1);
		if ((alpha3 > alpha2 && alpha3 > alpha1) || (alpha3 < alpha2 && alpha3 < alpha1)) {
			return beta;
		} else {
			return PI - beta;
		}
	}

	public static double normalizeAngle(double a) {
		a %= 2*PI;
		if (a > PI/2) {
			return a - PI;
		} else if (a < -PI/2) {
			return PI + a;
		} else {
			return a;
		}
	}

	/**
	 * Returns the angle between v1 and v2 in the range ]-pi/2, pi/2]. 
	 * Where the sign is the sign of the determinant |N v1 v2|. 
	 * @param v1
	 * @param v2
	 * @param N
	 * @return
	 */
	public static double getSignedAngle(double[] N, double[] v1, double[] v2) {
		double[][] T = {N, v1, v2};
		double sign = Math.signum(Rn.determinant(T));
		double alpha = Rn.euclideanAngle(v1, v2);
		if (alpha > PI/2) {
			alpha = -(PI - alpha);
		}
		return sign * alpha;
	}

	
	public static Map<Integer, Integer> createUndirectedEdgeMap(HalfEdgeDataStructure<?, ?, ?> hds) {
		Map<Integer, Integer> result = new HashMap<Integer, Integer>();
		Integer i = 0;
		for (Edge<?,?,?> e : hds.getPositiveEdges()) {
			result.put(e.getIndex(), i);
			result.put(e.getOppositeEdge().getIndex(), i);
			i++;
		}
		return result;
	}

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> double calculateAngleSumFromAlphas(V v, Map<E, Double> alphaMap) {
		double sum = 0.0;
		for (E e : HalfEdgeUtils.incomingEdges(v)) {
			sum += IsothermicUtility.getOppositeAlpha(e.getPreviousEdge(), alphaMap);
		}
		return sum;
	}

	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> double calculateAngleSumFromBetas(V v, Map<E, Double> betaMap) {
		double sum = 0.0;
		for (E e : HalfEdgeUtils.incomingEdges(v)) {
			if (e.getLeftFace() == null) {
				continue;
			}
			sum += betaMap.get(e.getPreviousEdge());
		}
		return sum;
	}
	
	
	public static int[] getPETScNonZeros(CoHDS hds, Functional<CoVertex, CoEdge, CoFace> fun){
		int [][] sparseStucture = fun.getNonZeroPattern(hds);
		int [] nnz = new int[sparseStucture.length];
		for(int i = 0; i < nnz.length; i++){
			nnz[i] = sparseStucture[i].length;
		}
		return nnz;
	}
	
	public static Map<CoFace, Double> calculateCirclePatternRhos(CoHDS hds, Map<CoEdge, Double> thetaMap, Map<CoFace, Double> phiMap) {
		// check pre-conditions
		System.out.println("Curvatures: ----------");
		double boundarySum = 0;
		for (CoVertex v : hds.getVertices()) {
			double Phi = 0.0;
			for (CoEdge e : HalfEdgeUtils.incomingEdges(v)) {
				Phi += thetaMap.get(e);
			}
			if (HalfEdgeUtils.isBoundaryVertex(v)) {
				boundarySum += 2*PI - Phi;
			}
			System.out.println(v + ": " + Phi/PI);
		}
		System.out.println("Boundary Sum/PI: " + boundarySum/PI);
		
		for (CoFace f : hds.getFaces()) {
			System.out.println(f + ": " + phiMap.get(f)/PI);
		}
		
		int dim = hds.numFaces();
		Vec rho = new Vec(dim);
		rho.zeroEntries();
		rho.assemble();
		
		CPEuclideanApplication app = new CPEuclideanApplication(hds, thetaMap, phiMap);
		app.setInitialSolutionVec(rho);
		Mat H = Mat.createSeqAIJ(dim, dim, PETSc.PETSC_DEFAULT, getPETScNonZeros(hds, app.getFunctional()));
		H.assemble();
		app.setHessianMat(H, H);
		
		Tao tao = new Tao(Method.NTR);
		tao.setApplication(app);
		tao.setMaximumIterates(200);
		tao.setTolerances(1E-15, 0, 0, 0);
		tao.setGradientTolerances(1E-15, 0, 0);
		tao.solve();
		System.out.println(tao.getSolutionStatus());
		
		Map<CoFace, Double> rhos = new HashMap<CoFace, Double>();
		for (CoFace f : hds.getFaces()) {
			rhos.put(f, rho.getValue(f.getIndex()));
		}
		return rhos;
	}


	public static Map<CoEdge, Double> calculateBetasFromAlphas(CoHDS hds, Map<CoEdge, Double> alphaMap) {
		Map<CoEdge, Double> betaMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : hds.getEdges()) {
			if (e.getLeftFace() != null) {
				double a1 = alphaMap.get(e);
				double a2 = alphaMap.get(e.getNextEdge());
				double a3 = alphaMap.get(e.getPreviousEdge());
				double beta = calculateTriangleAngle(a2, a3, a1);
				betaMap.put(e, beta);
			}
		}
		return betaMap;
	}
	
	public static Map<CoEdge, Double> calculateThetasFromBetas(CoHDS hds, Map<CoEdge, Double> betaMap) {
		Map<CoEdge, Double> thetas = new HashMap<CoEdge, Double>();
		for (CoEdge e : hds.getPositiveEdges()) {
			double theta = PI;
			if (e.getRightFace() != null) {
				double betaRight = betaMap.get(e.getOppositeEdge());
				theta -= betaRight;
			}
			if (e.getLeftFace() != null) {
				double betaLeft = betaMap.get(e);
				theta -= betaLeft;
			}
			thetas.put(e, theta);
			thetas.put(e.getOppositeEdge(), theta);
		}
		return thetas;
	}

	
	public static Map<CoFace, Double> calculatePhisFromBetas(CoHDS hds, Map<CoEdge, Double> betaMap) {
		Map<CoFace, Double> phiMap = new HashMap<CoFace, Double>();
		for (CoFace f : hds.getFaces()) {
			if (HalfEdgeUtils.isInteriorFace(f)) {
				phiMap.put(f, 2*PI);
			} else {
				double Phi = 2*PI;
				for (CoEdge e : HalfEdgeUtils.boundaryEdges(f)) {
					if (e.getRightFace() == null) {
						double beta = betaMap.get(e);
						Phi -= 2*beta;
					}
				}
				phiMap.put(f, Phi);
			}
		}
		return phiMap;
	}
	

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> double getOppositeAlpha(E e, Map<E, Double> alphaMap) {
		Double eA = alphaMap.get(e);
		if (eA == null) {
			eA = alphaMap.get(e.getOppositeEdge());
		}
		Double eNextA = alphaMap.get(e.getNextEdge());
		if (eNextA == null) {
			eNextA = alphaMap.get(e.getOppositeEdge().getNextEdge());
		}
		Double ePrevA = alphaMap.get(e.getPreviousEdge());
		if (ePrevA == null) {
			ePrevA = alphaMap.get(e.getOppositeEdge().getPreviousEdge());
		}
		return calculateTriangleAngle(eNextA, ePrevA, eA);
	}
	
}
