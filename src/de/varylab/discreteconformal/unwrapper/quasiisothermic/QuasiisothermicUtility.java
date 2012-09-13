package de.varylab.discreteconformal.unwrapper.quasiisothermic;

import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryVertex;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.atan;
import static java.lang.Math.log;
import static java.lang.Math.sin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
import de.jtem.halfedgetools.adapter.Adapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.CurvatureFieldMin;
import de.jtem.halfedgetools.adapter.type.Normal;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.jtem.halfedgetools.adapter.type.generic.BaryCenter4d;
import de.jtem.halfedgetools.adapter.type.generic.EdgeVector;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.halfedgetools.functional.Functional;
import de.jtem.halfedgetools.plugin.algorithm.vectorfield.EdgeVectorFieldMaxAdapter;
import de.jtem.jpetsc.InsertMode;
import de.jtem.jpetsc.KSP;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.MatStructure;
import de.jtem.jpetsc.PETSc;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.types.OppositeAngle;
import de.varylab.discreteconformal.unwrapper.ConesUtility;

public class QuasiisothermicUtility {

	static {
		Tao.Initialize();
	}
	
	@OppositeAngle
	public static class OppositeAnglesAdapter <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> extends AbstractAdapter<Double> {
		
		private Map<E, Double>
			angleMap = new HashMap<E, Double>();
		
		public OppositeAnglesAdapter(Map<E, Double> betaMap) {
			super(Double.class, true, true);
			this.angleMap = betaMap;
		}
		
		@Override
		public <
			VV extends Vertex<VV, EE, FF>,
			EE extends Edge<VV, EE, FF>,
			FF extends Face<VV, EE, FF>
		> Double getE(EE e, AdapterSet a) {
			if (!angleMap.containsKey(e)) {
				return 0.0;
			}
			return angleMap.get(e);
		}
		
		@SuppressWarnings("unchecked")
		@Override
		public <
			VV extends Vertex<VV, EE, FF>,
			EE extends Edge<VV, EE, FF>,
			FF extends Face<VV, EE, FF>
		> void setE(EE e, Double value, AdapterSet a) {
			angleMap.put((E)e, value);
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
		alpha1 = QuasiisothermicUtility.normalizeAngle(alpha1);
		alpha2 = QuasiisothermicUtility.normalizeAngle(alpha2);
		alpha3 = QuasiisothermicUtility.normalizeAngle(alpha3);
		double beta = abs(alpha2 - alpha1);
		if ((alpha3 > alpha2 && alpha3 > alpha1) || (alpha3 < alpha2 && alpha3 < alpha1)) {
			return beta;
		} else {
			return PI - beta;
		}
	}
	

	public static double normalizeAngle(double a) {
		a %= PI;
		if (a > PI/2) {
			return a - PI;
		} else if (a < -PI/2) {
			return PI + a;
		} else {
			return a;
		}
	}

	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void subdivideFaceSingularities(HDS hds, Map<E, Double> alpha, AdapterSet a) {
		Set<F> faces = new HashSet<F>(hds.getFaces());
		for (F f : faces) {
			double index = getSingularityIndexF(f, alpha, a);
			if (Math.abs(index) < 0.25) continue;
			
			// calulate alpha rotations for edges
			Map<E, Double> alphaRot = new HashMap<E, Double>();
			for (E e : HalfEdgeUtils.boundaryEdges(f)) {
				double ar = alphaRotationE(e, alpha, a);
				alphaRot.put(e, ar);
			}
			// subdivide and interpolate according to the transport
			double[] p = a.getD(BaryCenter4d.class, f);
			System.out.print("subdividing face " + f);
			V v = TopologyAlgorithms.splitFace(f);
			System.out.println("with index " + index + " -> " + v);
			a.set(Position.class, v, p);
			for (E e : HalfEdgeUtils.incomingEdges(v)) {
				E ep = e.getPreviousEdge();
				double ar = alphaRot.get(ep);
				double a1 = alpha.get(e.getPreviousEdge());
				double newAlpha = a1 + ar/2;
				alpha.put(e, newAlpha);
				alpha.put(e.getOppositeEdge(), newAlpha);
			}
		}
	}
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double getSingularityIndexF(F f, Map<E, Double> alpha, AdapterSet a) {
		return Math.round(2 * alphaRotationF(f, alpha, a) / (2.0 * Math.PI)) / 2.0;
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double getSingularityIndexV(V v, Map<E, Double> alpha, AdapterSet a) {
		return Math.round(2 * alphaRotationV(v, alpha, a) / (2.0 * Math.PI)) / 2.0;
	}
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> int getSingularityIndexFromAdapterSet(F f, AdapterSet a) {
		return (int)(Math.round(2 * alphaRotationFromAdapterSet(f, a) / (2.0 * Math.PI)) / 2.0);
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

	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<E, Integer> createSolverEdgeIndexMap(HDS hds, boolean excludeBoundary) {
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
	
	public static  <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<V, Integer> createSolverVertexIndexMap(HDS hds) {
		Map<V, Integer> result = new HashMap<V, Integer>();
		Integer i = 0;
		for (V v : hds.getVertices()) {
			if (isBoundaryVertex(v)) { 
				result.put(v, -1);
			} else {
				result.put(v, i);
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
			sum += QuasiisothermicUtility.getOppositeBeta(e.getPreviousEdge(), alphaMap);
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
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> int[] getPETScNonZeros(HDS hds, Functional<V, E, F> fun){
		int [][] sparseStucture = fun.getNonZeroPattern(hds);
		int [] nnz = new int[sparseStucture.length];
		for(int i = 0; i < nnz.length; i++){
			nnz[i] = sparseStucture[i].length;
		}
		return nnz;
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<E, Double> calculateBetasFromAlphas(HDS hds, Map<E, Double> alphaMap) {
		Map<E, Double> betaMap = new HashMap<E, Double>();
		for (E e : hds.getEdges()) {
			if (e.getLeftFace() != null) {
				double a1 = alphaMap.get(e);
				double a2 = alphaMap.get(e.getNextEdge());
				double a3 = alphaMap.get(e.getPreviousEdge());
				double beta = angleOrientation(a2,a3,a1)*calculateBeta(a2, a3, a1);
				betaMap.put(e, beta);
			}
		}
		return betaMap;
	}
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void checkTriangleAngles(HDS hds, Map<E, Double> betaMap) {
		double EPS = 1E-5;
		for (F f : hds.getFaces()) {
			double sum = 0.0;
			for (E e : HalfEdgeUtils.boundaryEdges(f)) {
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
		for (V v : hds.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) {
				continue;
			}
			double sum = 0.0;
			for (E e : HalfEdgeUtils.incomingEdges(v)) {
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
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void createDelaunayAngleSystem(HDS hds, Map<E, Double> betaMap) {
		OppositeAnglesAdapter<V, E, F> anglesAdapter = new OppositeAnglesAdapter<V, E, F>(betaMap);
		AdapterSet a = new AdapterSet(anglesAdapter);
		QuasiisothermicDelaunay.constructDelaunay(hds, a);
	}
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<E, Double> calculateThetasFromBetas(HDS hds, Map<E, Double> betaMap) {
		Map<E, Double> thetas = new HashMap<E, Double>();
		for (E e : hds.getPositiveEdges()) {
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

	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<F, Double> calculatePhisFromBetas(HDS hds, Map<E, Double> betaMap) {
		Map<F, Double> phiMap = new HashMap<F, Double>();
		for (F f : hds.getFaces()) {
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
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void layoutEars(HDS hds, Map<E, Double> alphaMap) {
//		List<E> earEdges = findEarsEdge(hds);
//		for (E ear : earEdges) {
//			E base = ear.getPreviousEdge();
//			V v = ear.getTargetVertex();
//			double[] p2 = base.getTargetVertex().T;
//			double[] p1 = base.getStartVertex().T;
//			double lBase = Pn.distanceBetween(p2, p1, Pn.EUCLIDEAN);
//			double lNew = getEdgeLength(ear, lBase, alphaMap);
//			
//		}
	}
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void alignLayout(HDS hds, Map<E, Double> alphaMap, AdapterSet a) {
		 List<E> b = HalfEdgeUtils.boundaryEdges(hds);
		 List<E> earEdges = findEarsEdge(hds);
		 for (E e : earEdges) {
			 b.remove(e);
			 b.remove(e.getNextEdge());
		 }
		 E refEdge = hds.getEdge(0);
		 if (b.size() != 0) {
			 refEdge = b.get(0);
		 }
		 double alpha = alphaMap.get(refEdge);
		 double[] s = a.getD(TexturePosition4d.class, refEdge.getStartVertex());
		 double[] t = a.getD(TexturePosition4d.class, refEdge.getTargetVertex());
		 double angle = atan((s[1] - t[1]) / (s[0] - t[0]));
		 double difAngle = alpha - angle;
		 MatrixBuilder mb = MatrixBuilder.euclidean();
		 mb.rotate(difAngle, new double[] {0,0,1});
		 Matrix T = mb.getMatrix();
		 for (V v : hds.getVertices()) {
			 double[] vt = a.getD(TexturePosition4d.class, v);
			 T.transformVector(vt);
			 a.set(TexturePosition.class, v, vt);
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
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void cutConesToBoundary(HDS hds, Map<E, Double> betaMap) {
		Set<V> innerVerts = new HashSet<V>(hds.getVertices());
		innerVerts.removeAll(HalfEdgeUtils.boundaryVertices(hds));
		for (V v : innerVerts) {
			if (!(v instanceof CoVertex)) throw new IllegalArgumentException("cutConesToBoundary() only works with CoVertex");
			CoVertex cv = (CoVertex)v;
			double sum = calculateAngleSumFromBetas(v, betaMap);
			if (abs(abs(sum) - 2*PI) > Math.PI/4) {
				System.out.println("angle sum " + sum);
				int index = (int)Math.round(sum / PI);
				cv.setTheta(index * PI);
				if((index % 2) != 0) {
					System.out.println("singularity: " + v + ", " + index + " pi");
				}
			} else {
				cv.setTheta(2 * PI);
			}
		}
		ConesUtility.cutMesh((CoHDS)hds);
	}

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<E, Double> calculateAlphasFromCurvature(AdapterSet a, HDS hds) {
		Map<E, Double> alphaMap = new HashMap<E, Double>();
		for(E e : hds.getEdges()) {
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

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double calculateAlpha(E e, AdapterSet a) {
		double[] N = a.getD(Normal.class, e);
		double[] Kmin = a.getD(CurvatureFieldMin.class, e);
		double[] E = a.getD(EdgeVector.class, e);
		return QuasiisothermicUtility.getSignedAngle(N, Kmin, E);
	}
	

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<F, Double> calculateOrientationFromAlphas(HDS hds, Map<E, Double> alphaMap) {
		Map<F, Double> orientationMap = new HashMap<F, Double>(hds.numFaces());
		for(F f : hds.getFaces()) {
			E e = f.getBoundaryEdge();
			double orientation = angleOrientation(
				alphaMap.get(e), 
				alphaMap.get(e.getNextEdge()), 
				alphaMap.get(e.getPreviousEdge())
			);
			orientationMap.put(f, orientation);
		}
		return orientationMap;
	}


	public static double angleOrientation(double alpha1, double alpha2,	double alpha3) {
		return (
			((alpha1 < alpha2) && (alpha2 < alpha3)) ||
			((alpha2 < alpha3) && (alpha3 < alpha1)) ||
			((alpha3 < alpha1) && (alpha1 < alpha2)) 
				)? -1.0 : 1.0;
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double alphaOrientationFromAdapterSet(F f, AdapterSet a) {
		E e = f.getBoundaryEdge();
		double alpha1 = calculateAlpha(e,a);
		double alpha2 = calculateAlpha(e.getNextEdge(),a);
		double alpha3 = calculateAlpha(e.getPreviousEdge(),a);
		return angleOrientation(alpha1, alpha2, alpha3);
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double alphaRotationFromAdapterSet(F f, AdapterSet a) {
		E se = f.getBoundaryEdge();
		E e = se;
		double rotation = 0.0;
		do {
			rotation += alphaRotationFromAdapterSet(e,a);
			e = e.getNextEdge();
		} while(e != se); 
		return rotation;
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double alphaRotationFromAdapterSet(V v, AdapterSet a) {
		E se = v.getIncomingEdge();
		E e = se;
		double rotation = 0.0;
		do {
			rotation += alphaRotationFromAdapterSet(e,a);
			e = e.getNextEdge().getOppositeEdge();
		} while(e != se);
		return rotation;
	}
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double alphaRotationFromAdapterSet(E e, AdapterSet a) {
		double alpha1 = calculateAlpha(e,a);
		double[] edge1 = a.getD(EdgeVector.class,e);
		double alpha2 = calculateAlpha(e.getNextEdge(),a);
		double[] edge2 = a.getD(EdgeVector.class,e.getNextEdge());
		double gamma = Rn.euclideanAngle(edge1, edge2);
		double rotation = normalizeAngle(gamma + alpha1 -alpha2 - Math.PI);
		return rotation;
	}
	
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double alphaRotationF(F f, Map<E, Double> alpha, AdapterSet a) {
		E se = f.getBoundaryEdge();
		E e = se;
		double rotation = 0.0;
		do {
			rotation += alphaRotationE(e, alpha, a);
			e = e.getNextEdge();
		} while(e != se); 
		return rotation;
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double alphaRotationV(V v, Map<E, Double> alpha, AdapterSet a) {
		E se = v.getIncomingEdge();
		E e = se;
		double rotation = 0.0;
		do {
			rotation += alphaRotationE(e, alpha, a);
			e = e.getNextEdge().getOppositeEdge();
		} while(e != se);
		return rotation;
	}
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double alphaRotationE(E e, Map<E, Double> alpha, AdapterSet a) {
		double alpha1 = alpha.get(e);
		double[] edge1 = a.getD(EdgeVector.class,e);
		double alpha2 = alpha.get(e.getNextEdge());
		double[] edge2 = a.getD(EdgeVector.class,e.getNextEdge());
		double gamma = Rn.euclideanAngle(edge1, edge2);
		double rotation = normalizeAngle(gamma + alpha1 -alpha2 - Math.PI);
		return rotation;
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<V, Double> calculateQuasiconformalFactors(HDS hds, Map<E, Double> edgeLengths, Map<E, Double> texLengths) {
		Map<E, Integer> eIndexMap = createSolverEdgeIndexMap(hds, false);
		int m = hds.numEdges() / 2;
		Vec u = new Vec(m);
		Vec l = new Vec(m);
		Mat A = Mat.createSeqAIJ(m, m, 2, null);
		for (E e : hds.getPositiveEdges()) {
			int eIndex = eIndexMap.get(e);
			int i = e.getStartVertex().getIndex();
			int j = e.getTargetVertex().getIndex();
			double l1 = 2*log(edgeLengths.get(e));
			double l2 = 2*log(texLengths.get(e));
			l.setValue(eIndex, l2 - l1, InsertMode.INSERT_VALUES);
			A.setValue(eIndex, i, 1.0, InsertMode.INSERT_VALUES);
			A.setValue(eIndex, j, 1.0, InsertMode.INSERT_VALUES);
		}
		A.assemble();
		KSP ksp = KSP.create();
		ksp.setOptionsPrefix("qcf_");
		PETSc.optionsSetValue("-qcf_ksp_type", "lsqr");

		ksp.setFromOptions();
//		ksp.setTolerances(1E-8, PETSc.PETSC_DEFAULT, PETSc.PETSC_DEFAULT, 25);
		ksp.setOperators(A, A, MatStructure.SAME_NONZERO_PATTERN);
		ksp.solve(l, u);
		System.out.println("ksp residual: " + ksp.getResidualNorm());
		Map<V, Double> result = new HashMap<V, Double>();
		for (V v : hds.getVertices()) {
			double ui = u.getValue(v.getIndex());
			result.put(v, ui);
		}
		return result;
	}


	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Adapter<double[]> createAlphaField(HDS hds, Map<E, Double> alpha, AdapterSet a, String name) {
		Map<E, double[]> vMap = new HashMap<E, double[]>();
		for (E e : hds.getPositiveEdges()) {
			if (!alpha.containsKey(e)) continue;
			double al = alpha.get(e);
			double[] ev = a.getD(EdgeVector.class, e);
			Rn.normalize(ev, ev);
			double[] n = a.getD(Normal.class, e);
			double[] vp = Rn.crossProduct(null, ev, n);
			Rn.times(ev, Math.cos(al), ev);
			Rn.times(vp, Math.sin(al), vp);
			double[] vec = Rn.add(null, ev, vp);
			vMap.put(e, vec);
			vMap.put(e.getOppositeEdge(), vec);
		}
		EdgeVectorFieldMaxAdapter aMax = new EdgeVectorFieldMaxAdapter(vMap, name);
		return aMax;
	}
	
}
