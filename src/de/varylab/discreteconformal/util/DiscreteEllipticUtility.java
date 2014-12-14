package de.varylab.discreteconformal.util;

import static java.lang.Math.abs;
import static java.lang.Math.signum;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.logging.Logger;

import com.wolfram.jlink.KernelLink;
import com.wolfram.jlink.MathLinkException;
import com.wolfram.jlink.MathLinkFactory;

import de.jreality.math.Pn;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.Unwrapper;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;


public class DiscreteEllipticUtility {

	private static Logger
		log = Logger.getLogger(DiscreteEllipticUtility.class.getName());
	
	/**
	 * Calculates the half-period ratio of an elliptic function. 
	 * The parameter domain is given by the texture coordinates and the identification
	 * defined in cutInfo.
	 * @param cutInfo The period lattice identification for the mesh 
	 * @return the half-period ratio tau
	 */
	public static Complex calculateHalfPeriodRatio(CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		Set<CoVertex> branchSet = cutInfo.getCopies(cutInfo.cutRoot);
		if (branchSet.size() != 4) {
			log.warning("domain torus has more than one vertex");
			return new Complex(0, 0);
		}
		Iterator<CoVertex> branchIt = branchSet.iterator();
		CoVertex v0 = branchIt.next();
		CoVertex v1 = branchIt.next();
		CoVertex v2 = branchIt.next();
		double[] t0 = v0.T;
		double[] t1 = v1.T;
		double[] t2 = v2.T;
		Complex z0 = new Complex(t0[0] / t0[3], t0[1] / t0[3]);
		Complex z1 = new Complex(t1[0] / t1[3], t1[1] / t1[3]);
		Complex z2 = new Complex(t2[0] / t2[3], t2[1] / t2[3]);
		Complex w1 = z1.minus(z0);
		Complex w2 = z2.minus(z0);
		Complex tau = w2.divide(w1);
		Complex tauNorm = normalizeModulus(tau);
		return tauNorm;
	}
	
	/**
	 * Move tau to the fundamental domain of flat tori
	 * @param tau
	 * @return
	 */
	public static Complex normalizeModulus(Complex t) {
		Complex tau = new Complex(t);
		int maxIter = 100;
		// move tau into its fundamental domain
		while ((abs(tau.re) > 0.5 || tau.im < 0 || tau.re < 0 || tau.abs() < 1) && --maxIter > 0) {
			if (abs(tau.re) > 0.5) {
				tau.re -= signum(tau.re);
			}
			if (tau.im < 0) {
				tau.im *= -1;
			}
			if (tau.re < 0) {
				tau.re *= -1;
			}
			if (tau.abs() < 1) {
				tau = tau.invert();
			}
		}
		return tau;
	}

	/**
	 * Calculates the half-period ratio for a given doubly covered branched sphere.
	 * First the elliptic function is approximated by calculating a conformally flat
	 * metric. Then the corresponding triangle layout in the complex plane is calculated.
	 * From this we obtain the periods by the coordinates of the indentified vertices
	 * along the cut paths. 
	 * @param hds A triangulated torus with vertices on the sphere.
	 * @return The half-period ratio of the elliptic function
	 */
	public static Complex calculateHalfPeriodRatio(CoHDS hds, CoVertex cutRoot, double tol, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) throws Exception {
		Unwrapper unwrapper = new EuclideanUnwrapperPETSc();
		unwrapper.setGradientTolerance(tol);
		unwrapper.setMaxIterations(1000);
		AdapterSet a = new ConformalAdapterSet();
		unwrapper.setCutRoot(cutRoot);
		unwrapper.unwrap(hds, 1, a);
		CuttingInfo<CoVertex, CoEdge, CoFace> info = unwrapper.getCutInfo();
		if (cutInfo != null) {
			cutInfo.cutRoot = info.cutRoot;
			cutInfo.edgeCutMap = info.edgeCutMap;
			cutInfo.vertexCopyMap = info.vertexCopyMap;
		}
		return calculateHalfPeriodRatio(info);
	}

	public static Map<CoVertex, CoVertex> generateEllipticImage(CoHDS hds, int numExtraPoints, Set<CoEdge> glueEdges, int... branchVertices) {
		return generateEllipticImage(hds, numExtraPoints, true, true, glueEdges, branchVertices);
	}
	

	/**
	 * Generate a doubly covered sphere with four branch points. This is combinatorially 
	 * a torus and can be interpreted as the image of an elliptic function.
	 * The branch points are picked from the given vertices, either the first 4 or as defined in branchVertices.
	 * The second sheet of the double cover is a copy of the first one glued along paths between each two branch points.
	 * @param hds The initial vertices of the image, must have at least four vertices. The edges and faces are ignored 
	 * @param numExtraPoints Number of extra random points added to the initial points
	 * @param useConvexHull use edges of convex hull or 
	 * @param glueEdges The edges 
	 * @param branchVertices
	 * @param returns the vertex-wise hyperelliptic involution. At branch points one of the vertices is removed.
	 */
	public static Map<CoVertex, CoVertex> generateEllipticImage(CoHDS hds, int numExtraPoints, boolean projectToSphere, boolean useConvexHull, Set<CoEdge> glueEdges, int... branchVertices) {
		if (numExtraPoints > 0 && !useConvexHull) {
			throw new IllegalArgumentException("cannot add extra points if useConvexHull is false");
		}
		if (hds.numVertices() < 4) {
			throw new RuntimeException("No branch point set in generateEllipticCurve()");
		} 
		if (branchVertices.length < 4) {
			branchVertices = new int[]{0, 1, 2, 3};
		}
		CoVertex v1 = hds.getVertex(branchVertices[0]);
		CoVertex v2 = hds.getVertex(branchVertices[1]);
		CoVertex v3 = hds.getVertex(branchVertices[2]);
		CoVertex v4 = hds.getVertex(branchVertices[3]);
	
		// additional points
		Random rnd = new Random();
		for (int i = 0; i < numExtraPoints; i++) {
			CoVertex v = hds.addNewVertex();
			v.P = new double[] {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian(), 1.0};
		}
		
		// project to sphere
		if (projectToSphere) {
			for (CoVertex v : hds.getVertices()) {
				Pn.setToLength(v.P, v.P, 1.0, Pn.EUCLIDEAN);
			}
		}
		
		// create new triangulation via convex hull
		if (useConvexHull) {
			ConformalAdapterSet a = new ConformalAdapterSet();
			for (CoEdge e : new HashSet<CoEdge>(hds.getEdges())) {
				hds.removeEdge(e);
			}
			for (CoFace f : new HashSet<CoFace>(hds.getFaces())) {
				hds.removeFace(f);
			}
			ConvexHull.convexHull(hds, a, 1E-8);
		}
		
		int vOffset = hds.numVertices();
		int eOffset = hds.numEdges();
		HalfEdgeUtils.copy(hds, hds);
		Map<CoVertex, CoVertex> involution = new HashMap<>();
		for (int i = 0; i < vOffset; i++) {
			CoVertex v = hds.getVertex(i);
			CoVertex vc = hds.getVertex(vOffset + i); 
			involution.put(v, vc);
			involution.put(vc, v);
			double[] p = v.P;
			vc.P = p.clone();
		}
		
		Set<CoVertex> path2Ends = new HashSet<CoVertex>();
		path2Ends.add(v3);
		path2Ends.add(v4);
		List<CoEdge> path1 = Search.bFS(v1, v2, path2Ends);
		Set<CoVertex> path1Vertices = PathUtility.getVerticesOnPath(path1);
		List<CoEdge> path2 = Search.bFS(v3, v4, path1Vertices);
		
		List<CoEdge> path1c = new LinkedList<CoEdge>();
		List<CoEdge> path2c = new LinkedList<CoEdge>();
		for (CoEdge e : path1) {
			path1c.add(hds.getEdge(eOffset + e.getIndex()));
		}
		for (CoEdge e : path2) {
			path2c.add(hds.getEdge(eOffset + e.getIndex()));
		}
		
		SurgeryUtility.cutAndGluePaths(path1, path1c);
		SurgeryUtility.cutAndGluePaths(path2, path2c);
		
		if (glueEdges != null) {
			glueEdges.clear();
			glueEdges.addAll(path1);
			glueEdges.addAll(path2);
			glueEdges.addAll(path1c);
			glueEdges.addAll(path2c);
		}
		return involution;
	}
	
	public static Complex calculateHalfPeriodRatioMathLink(double[] p1, double[] p2, double[] p3, double[] p4, KernelLink l) throws MathLinkException {
		// to the sphere
		Pn.setToLength(p1, p1, 1.0, Pn.EUCLIDEAN);
		Pn.setToLength(p2, p2, 1.0, Pn.EUCLIDEAN);
		Pn.setToLength(p3, p3, 1.0, Pn.EUCLIDEAN);
		Pn.setToLength(p4, p4, 1.0, Pn.EUCLIDEAN);
		Pn.dehomogenize(p1, p1);
		Pn.dehomogenize(p2, p2);
		Pn.dehomogenize(p3, p3);
		Pn.dehomogenize(p4, p4);
		// project stereographically
		Complex z1 = new Complex(p1[0] / (1 - p1[2]), p1[1] / (1 - p1[2]));
		Complex z2 = new Complex(p2[0] / (1 - p2[2]), p2[1] / (1 - p2[2]));
		Complex z3 = new Complex(p3[0] / (1 - p3[2]), p3[1] / (1 - p3[2]));
		Complex z4 = new Complex(p4[0] / (1 - p4[2]), p4[1] / (1 - p4[2]));
		Complex e1 = toInfinitZeroSum(z1, z1, z2, z3, z4);
		Complex e2 = toInfinitZeroSum(z2, z1, z2, z3, z4);
		Complex e3 = toInfinitZeroSum(z3, z1, z2, z3, z4);
		Complex g2 = e1.times(e2).plus(e2.times(e3)).plus(e3.times(e1)).times(-4);
		Complex g3 = e1.times(e2).times(e3).times(4);
		Complex[] w1w3 = invokeWeierstrassHalfperiods(g2, g3, l);
		Complex tau = w1w3[1].divide(w1w3[0]);
		return normalizeModulus(tau);
	}
	

	private static Complex toInfinitZeroSum(Complex z, Complex z1, Complex z2, Complex z3, Complex z4) {
		Complex t1 = z.minus(z4).invert();
		Complex t2 = z1.times(z2.plus(z3).minus(z4.times(2))).plus(z2.times(z3.minus(z4.times(2)))).plus(z4.times(z4.times(3).minus(z3.times(2))));
		Complex t3 = z4.minus(z1).times(z2.minus(z4)).times(z4.minus(z3)).times(3);
		return t1.minus(t2.divide(t3));
	}
	
	
	
	public static Complex[] invokeWeierstrassHalfperiods(Complex w2, Complex w3, KernelLink l) throws MathLinkException {
		l.setComplexClass(MLComplex.class);
		MLComplex[] g2g3 = new MLComplex[] {new MLComplex(w2), new MLComplex(w3)};
		l.putFunction("EvaluatePacket", 1);
		l.putFunction("WeierstrassHalfPeriods", 1);
		l.put(g2g3);
		l.endPacket();
		l.waitForAnswer();
		MLComplex[] r = (MLComplex[])l.getComplexArray1();
		return r;
	}
	
	
	
	public static void main(String[] args) throws Exception {
		NativePathUtility.set("native");
		String[] mlargs = new String[] {
			"-linkmode", "launch", 
			"-linkname", "\"C:\\Program Files\\Wolfram Research\\Mathematica\\7.0\\MathKernel.exe\" " + 
			"-mathlink"
		};
		KernelLink link = MathLinkFactory.createKernelLink(mlargs);
		link.discardAnswer();
		
		double[] p1 = {0.5257311122273196, 0.8506508082851774, -2.259622999145302E-9, 1.0};
		double[] p2 = {-7.865783539640066E-9, -0.525731111594899, -0.8506508086760348, 1.0};
		double[] p3 = {0.5257311164974198, -0.8506508056461103, -1.2759862959410077E-10, 1.0};
		double[] p4 = {5.046344158026151E-10, 0.5257311116954958, 0.8506508086138624, 1.0};
		Complex tau1 = calculateHalfPeriodRatioMathLink(p1, p2, p3, p4, link);
		System.out.println("tau1 = " + tau1);
		
		
		Complex g2 = new Complex(-0.0143277, 0.0158116);
		Complex g3 = new Complex(-0.000768009, 2.78129E-6);
		Complex[] w1w3 = invokeWeierstrassHalfperiods(g2, g3, link);
		Complex tau2 = w1w3[1].divide(w1w3[0]);
		System.out.println("tau2 = " + tau2);
	}
	
	
	
	public static class MLComplex extends Complex {

		private static final long 
			serialVersionUID = 1L;

		public MLComplex() {
			super();
		}

		public MLComplex(Complex u) {
			super(u);
		}

		public MLComplex(de.jtem.mfc.field.Field.Complex u) {
			super(u);
		}

		public MLComplex(double aReal, double aImag) {
			super(aReal, aImag);
		}

		public MLComplex(double aReal) {
			super(aReal);
		}

		public double re() {
			return re;
		}
		
		public double im() {
			return im;
		}
		
	}
	

	
}
