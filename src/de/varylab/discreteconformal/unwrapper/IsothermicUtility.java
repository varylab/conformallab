package de.varylab.discreteconformal.unwrapper;

import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;
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
import de.varylab.discreteconformal.functional.CPEuclideanFunctional.Phi;
import de.varylab.discreteconformal.functional.CPEuclideanFunctional.Theta;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.CPLayoutAdapters.Radius;
import de.varylab.discreteconformal.unwrapper.CPLayoutAdapters.Rho;
import de.varylab.discreteconformal.unwrapper.CPLayoutAdapters.XYFace;
import de.varylab.discreteconformal.unwrapper.CPLayoutAdapters.XYVertex;
import de.varylab.discreteconformal.unwrapper.isothermic.IsothermicLayout;
import de.varylab.discreteconformal.unwrapper.numerics.CPEuclideanApplication;

public class IsothermicUtility {

	public static class CPFunctionalAdapters implements Theta<CoEdge>, Phi<CoFace> {
	
		public Map<CoEdge, Double>
			thetaMap = new HashMap<CoEdge, Double>();
		public Map<CoFace, Double>
			phiMap = new HashMap<CoFace, Double>();
		
		@Override
		public double getPhi(CoFace f) {
			if (HalfEdgeUtils.isInteriorFace(f)) {
				return 2 * PI;
			} else {
				double Phi = 2*PI;
				for (CoEdge e : HalfEdgeUtils.boundaryEdges(f)) {
					if (!HalfEdgeUtils.isBoundaryEdge(e)) {
						continue;
					}
					double eStar = PI - getTheta(e);
					Phi -= 2*eStar;
				}
				return Phi;
			}
		}
	
		@Override
		public void setPhi(CoFace f, double phi) {
			phiMap.put(f, phi);
		}
	
		@Override
		public double getTheta(CoEdge e) {
			if (!thetaMap.containsKey(e)) {
				return PI;
			}
			return thetaMap.get(e);
		}
	
		@Override
		public void setTheta(CoEdge e, double theta) {
			thetaMap.put(e, theta);
		}
		
	}

	public static class CPLayoutAdapters implements XYVertex<CoVertex>, XYFace<CoFace>, Rho<CoFace>, Radius<CoFace> {
	
		public Vec rho = null;
		
		public CPLayoutAdapters(Vec rho) {
			super();
			this.rho = rho;
		}
	
		@Override
		public double getRadius(CoFace v) {
			return 0;
		}
	
		@Override
		public void setRadius(CoFace v, double r) {
		}
	
		@Override
		public double getRho(CoFace f) {
			return rho.getValue(f.getIndex());
		}
	
		@Override
		public void setRho(CoFace f, double rho) {
			this.rho.setValue(f.getIndex(), rho, INSERT_VALUES);
		}
	
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
			sum += IsothermicLayout.getOppositeAlpha(e.getPreviousEdge(), alphaMap);
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
	
	public static Vec calculateCirclePatternRadii(CoHDS hds, Map<CoEdge, Double> alphaMap, IsothermicUtility.CPFunctionalAdapters funAdapters) {
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
			funAdapters.setTheta(e, theta);
			funAdapters.setTheta(e.getOppositeEdge(), theta);
		}
		
		// check pre-conditions
		double boundarySum = 0;
		for (CoVertex v : hds.getVertices()) {
			double Phi = calculateAngleSumFromBetas(v, betaMap);
			if (HalfEdgeUtils.isBoundaryVertex(v)) {
				boundarySum += PI - Phi;
			} else {
				System.out.println(v + ": " + Phi);
			}
		}
		System.out.println("Boundary Sum/PI: " + boundarySum/PI);
		
		int dim = hds.numFaces();
		Vec rho = new Vec(dim);
		rho.zeroEntries();
		rho.assemble();
		
		CPEuclideanApplication app = new CPEuclideanApplication(hds, funAdapters, funAdapters);
		app.setInitialSolutionVec(rho);
		Mat H = Mat.createSeqAIJ(dim, dim, PETSc.PETSC_DEFAULT, getPETScNonZeros(hds, app.getFunctional()));
		H.assemble();
		app.setHessianMat(H, H);
		
		Tao tao = new Tao(Method.NTR);
		tao.setFromOptions();
		tao.setApplication(app);
		tao.setMaximumIterates(20);
		tao.setTolerances(1E-15, 0, 0, 0);
		tao.setGradientTolerances(1E-15, 0, 0);
		tao.solve();
		System.out.println(tao.getSolutionStatus());
		
		return app.getSolutionVec();
	}
	
}
