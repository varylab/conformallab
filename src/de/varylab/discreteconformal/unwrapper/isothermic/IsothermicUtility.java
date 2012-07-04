package de.varylab.discreteconformal.unwrapper.isothermic;

import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryVertex;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.atan;
import static java.lang.Math.sin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.vecmath.Point2d;

import de.jreality.math.Matrix;
import de.jreality.math.MatrixBuilder;
import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.CurvatureFieldMin;
import de.jtem.halfedgetools.adapter.type.Normal;
import de.jtem.halfedgetools.adapter.type.generic.EdgeVector;
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
import de.varylab.discreteconformal.heds.adapter.types.OppositeAngle;
import de.varylab.discreteconformal.unwrapper.ConesUtility;
import de.varylab.discreteconformal.unwrapper.circlepattern.CPLayoutAdapters.XYFace;
import de.varylab.discreteconformal.unwrapper.circlepattern.CPLayoutAdapters.XYVertex;
import de.varylab.discreteconformal.unwrapper.circlepattern.CPLayoutAlgorithm;
import de.varylab.discreteconformal.unwrapper.numerics.CPEuclideanApplication;

public class IsothermicUtility {

	public static class CPLayoutAdapters implements XYVertex<CoVertex>, XYFace<CoFace> {

		private Map<CoFace, double[]>
			centerMap = new HashMap<CoFace, double[]>();
		
		@Override
		public Point2d getXY(CoFace f, Point2d xy) {
			if (!centerMap.containsKey(f)) {
				centerMap.put(f, new double[2]);
			}
			return new Point2d(centerMap.get(f));
		}
	
		@Override
		public void setXY(CoFace f, Point2d xy) {
			centerMap.put(f, new double[] {xy.x, xy.y});
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
	
	@OppositeAngle
	public static class OppositeAnglesAdapter extends AbstractAdapter<Double> {
		
		private Map<CoEdge, Double>
			angleMap = new HashMap<CoEdge, Double>();
		
		public OppositeAnglesAdapter(Map<CoEdge, Double> betaMap) {
			super(Double.class, true, true);
			this.angleMap = betaMap;
		}
		
		@Override
		public <
			V extends Vertex<V, E, F>,
			E extends Edge<V, E, F>,
			F extends Face<V, E, F>
		> Double getE(E e, AdapterSet a) {
			if (!angleMap.containsKey(e)) {
				return 0.0;
			}
			return angleMap.get(e);
		}
		
		@Override
		public <
			V extends Vertex<V, E, F>,
			E extends Edge<V, E, F>,
			F extends Face<V, E, F>
		> void setE(E e, Double value, AdapterSet a) {
			angleMap.put((CoEdge)e, value);
		}
		
		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return CoEdge.class.isAssignableFrom(nodeClass);
		}
		
	}
	

	/**
	 * Calculate the angle between the edges that belong to alpha1 and alpha2.
	 * @param alpha1
	 * @param alpha2
	 * @param alpha3
	 * @return
	 */
	public static double calculateBeta(double alpha1, double alpha2, double alpha3) {
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

	
	public static  <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<E, Integer> createSolverIndexMap(HDS hds, boolean excludeBoundary) {
		Map<E, Integer> result = new HashMap<E, Integer>();
		Integer i = 0;
		for (E e : hds.getPositiveEdges()) {
			if (((isBoundaryVertex(e.getStartVertex()) || 
				isBoundaryVertex(e.getTargetVertex())) && excludeBoundary)
			) {
				result.put(e, -1);
				result.put(e.getOppositeEdge(), -1);
			} else {
				result.put(e, i);
				result.put(e.getOppositeEdge(), i);
				i++;
			}
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
			sum += IsothermicUtility.getOppositeBeta(e.getPreviousEdge(), alphaMap);
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
		double thetaStarSum = 0.0;
		for (CoEdge e : hds.getPositiveEdges()) {
			double thStar = PI-thetaMap.get(e);
			thetaStarSum += 2*thStar;
		}
		System.out.println("theta star sum: " + thetaStarSum/(2*PI));
		for (CoEdge e : hds.getPositiveEdges()) {
			double theta = thetaMap.get(e);
			if (theta < 0) {
				System.err.println("negative theta at " + e + ": " + theta);
			}
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
				double beta = calculateBeta(a2, a3, a1);
				betaMap.put(e, beta);
			}
		}
		return betaMap;
	}
	
	
	public static void checkTriangleAngles(CoHDS hds, Map<CoEdge, Double> betaMap) {
		double EPS = 1E-5;
		for (CoFace f : hds.getFaces()) {
			double sum = 0.0;
			for (CoEdge e : HalfEdgeUtils.boundaryEdges(f)) {
				double beta = betaMap.get(e);
				if (beta < 0) {
					throw new RuntimeException("Negative angle at edge " + e + ": " + beta);
				}
				sum += beta;
			}
			if (Math.abs(sum - PI) > EPS) {
				throw new RuntimeException("Angle sum at " + f + ": " + sum);
			}
		}
		for (CoVertex v : hds.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) {
				continue;
			}
			double sum = 0.0;
			for (CoEdge e : HalfEdgeUtils.incomingEdges(v)) {
				double beta = betaMap.get(e.getPreviousEdge());
				sum += beta;
			}
			if (Math.abs(sum - 2*PI) > EPS) {
				throw new RuntimeException("Angle sum at " + v + ": " + sum);
			}
		}
	}
	
	
	/**
	 * invokes the Delaunay flip algorithm on hds and modifies the angles in betaMap 
	 * accordingly. The local Delaunay condition is derived from the opposite angles
	 * of an edge.
	 * @param hds
	 * @param betaMap
	 */
	public static void createDelaunayAngleSystem(CoHDS hds, Map<CoEdge, Double> betaMap) {
		OppositeAnglesAdapter anglesAdapter = new OppositeAnglesAdapter(betaMap);
		AdapterSet a = new AdapterSet(anglesAdapter);
		IsothermicDelaunay.constructDelaunay(hds, a);
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
			phiMap.put(f, 2*PI);
		}
		return phiMap;
	}
	

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> double getOppositeBeta(E e, Map<E, Double> alphaMap) {
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
		return calculateBeta(eNextA, ePrevA, eA);
	}
	
	
	public static void layoutEars(CoHDS hds, Map<CoEdge, Double> alphaMap) {
//		List<CoEdge> earEdges = findEarsEdge(hds);
//		for (CoEdge ear : earEdges) {
//			CoEdge base = ear.getPreviousEdge();
//			CoVertex v = ear.getTargetVertex();
//			double[] p2 = base.getTargetVertex().T;
//			double[] p1 = base.getStartVertex().T;
//			double lBase = Pn.distanceBetween(p2, p1, Pn.EUCLIDEAN);
//			double lNew = getEdgeLength(ear, lBase, alphaMap);
//			
//		}
	}
	
	
	public static void alignLayout(CoHDS hds, Map<CoEdge, Double> alphaMap) {
		 List<CoEdge> b = HalfEdgeUtils.boundaryEdges(hds);
		 List<CoEdge> earEdges = findEarsEdge(hds);
		 for (CoEdge e : earEdges) {
			 b.remove(e);
			 b.remove(e.getNextEdge());
		 }
		 CoEdge refEdge = hds.getEdge(0);
		 if (b.size() != 0) {
			 refEdge = b.get(0);
		 }
		 double alpha = alphaMap.get(refEdge);
		 double[] s = refEdge.getStartVertex().T;
		 double[] t = refEdge.getTargetVertex().T;
		 double angle = atan((s[1] - t[1]) / (s[0] - t[0]));
		 double difAngle = alpha - angle;
		 MatrixBuilder mb = MatrixBuilder.euclidean();
		 mb.rotate(difAngle, new double[] {0,0,1});
		 Matrix T = mb.getMatrix();
		 for (CoVertex v : hds.getVertices()) {
			 T.transformVector(v.T);
		 }
	}
	
	

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> List<E> findEarsEdge(HDS hds){
		ArrayList<E> result = new ArrayList<E>();
		for (E e : hds.getEdges()){
			if (!HalfEdgeUtils.isInteriorEdge(e) && e.getLeftFace() == null) {
				if (e.getRightFace() == e.getNextEdge().getRightFace()) {
					result.add(e);
				}
			}
		}
		return result;
	}

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> double getEdgeLength(E e, double prevEdgeLength, Map<E, Double> alphaMap) {
		double ea = getOppositeBeta(e, alphaMap);
		double prevAngle = getOppositeBeta(e.getPreviousEdge(), alphaMap);
		return prevEdgeLength  * sin(ea) / sin(prevAngle);
	}

	
	/**
	 * Cuts from cone points (beta sum != 2PI) to the boundary of a mesh
	 * @param hds
	 * @param betaMap
	 */
	public static void cutConesToBoundary(CoHDS hds, Map<CoEdge, Double> betaMap) {
		Set<CoVertex> innerVerts = new HashSet<CoVertex>(hds.getVertices());
		innerVerts.removeAll(HalfEdgeUtils.boundaryVertices(hds));
		for (CoVertex v : innerVerts) {
			double sum = calculateAngleSumFromBetas(v, betaMap);
			if (abs(sum - 2*PI) > Math.PI/4) {
				int index = (int)Math.round(sum / PI);
				v.setTheta(index * PI);
				System.out.println("singularity: " + v + ", " + index + " pi");
			} else {
				v.setTheta(2 * PI);
			}
		}
		ConesUtility.cutMesh(hds);
	}

	public static Map<CoEdge, Double> calculateAlphasFromCurvature(AdapterSet a, CoHDS hds) {
		Map<CoEdge, Double> alphaMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : hds.getEdges()) {
			double[] N = a.getD(Normal.class, e);
			double[] Kmin = a.getD(CurvatureFieldMin.class, e);
			double[] E = a.getD(EdgeVector.class, e);
			try {
				double ae = getSignedAngle(N, Kmin, E);
				alphaMap.put(e, ae);
				alphaMap.put(e.getOppositeEdge(), ae);
			} catch (Exception e2) { 
				System.out.println("check");
			}
		}
		return alphaMap;
	}

	public static void doCirclePatternLayout(CoHDS hds, Map<CoEdge, Double> thetaMap, Map<CoFace, Double> rhoMap) {
		CPLayoutAdapters layoutAdapters = new CPLayoutAdapters();
		CPLayoutAlgorithm<CoVertex, CoEdge, CoFace>
			layout = new CPLayoutAlgorithm<CoVertex, CoEdge, CoFace>(
				layoutAdapters, layoutAdapters, rhoMap, thetaMap
			);
		layout.execute(hds);
	}
	
	
}
