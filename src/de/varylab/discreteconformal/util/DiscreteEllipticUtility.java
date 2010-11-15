package de.varylab.discreteconformal.util;

import static java.lang.Math.abs;
import static java.lang.Math.signum;
import geom3d.Point;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;

import no.uib.cipr.matrix.Vector;

import com.wolfram.jlink.KernelLink;
import com.wolfram.jlink.MathLinkException;
import com.wolfram.jlink.MathLinkFactory;

import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.TypedAdapterSet;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.unwrapper.EuclideanLayout;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.Unwrapper;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.Search.DefaultWeightAdapter;


public class DiscreteEllipticUtility {

	
	/**
	 * Calculates the half-period ratio of an elliptic function. 
	 * The parameter domain is given by the texture coordinates and the identification
	 * defined in cutInfo.
	 * @param cutInfo The period lattice identification for the mesh 
	 * @return the half-period ratio tau
	 */
	public static Complex calculateHalfPeriodRatio(CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		CoVertex v0 = cutInfo.cutRoot;
		CoVertex v1 = cutInfo.vertexCopyMap.get(v0);
		if (v1 == null) throw new RuntimeException("Connot calculate modulus. No cut-root copies found");
		CoVertex v2 = cutInfo.vertexCopyMap.get(v1);
		if (v2 == null) throw new RuntimeException("Connot calculate modulus. No second cut-root copy found");
		Point t0 = v0.getTextureCoord();
		Point t1 = v1.getTextureCoord();
		Point t2 = v2.getTextureCoord();
		Complex z0 = new Complex(t0.x() / t0.z(), t0.y() / t0.z());
		Complex z1 = new Complex(t1.x() / t1.z(), t1.y() / t1.z());
		Complex z2 = new Complex(t2.x() / t2.z(), t2.y() / t2.z());
		Complex w1 = z1.minus(z0);
		Complex w2 = z2.minus(z0);
		Complex tau = w2.divide(w1);
		
		int maxIter = 100;
		// move tau into its fundamental domain
		while ((abs(tau.re) > 0.5 || tau.im < 0 || tau.abs() < 1) && --maxIter > 0) {
			if (abs(tau.re) > 0.5) {
				tau.re -= signum(tau.re);
			}
			if (tau.im < 0) {
				tau = tau.times(-1);
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
	 * along the cutted paths. 
	 * @param hds A triangulated torus with vertices on the sphere.
	 * @return The half-period ratio of the elliptic function
	 */
	public static Complex calculateHalfPeriodRatio(CoHDS hds, double tol) {
		Unwrapper unwrapper = new EuclideanUnwrapperPETSc();
		unwrapper.setGradientTolerance(tol);
		unwrapper.setMaxIterations(500);
		Vector u = null;
		try {
			u = unwrapper.unwrap(hds, new AdapterSet());
		} catch (Exception e) {
			e.printStackTrace();
			return new Complex();
		}
		DefaultWeightAdapter<CoEdge> w = new DefaultWeightAdapter<CoEdge>();
		CoVertex cutRoot = hds.getVertex(0);
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = CuttingUtility.cutTorusToDisk(hds, cutRoot, w);
		EuclideanLayout.doLayout(hds, u);
		return calculateHalfPeriodRatio(cutInfo);
	}

	

	/**
	 * Generate a doubly covered sphere with four branch points. This is combinatorially 
	 * a torus and can be interpreted as the image of an elliptic function.
	 * The branch points are picked from the given vertices, either the first 4 or as defined in branchVertices.
	 * The second sheet of the double cover is a copy of the first one glued along paths between each two branch points.
	 * @param hds The initial vertices of the image, must have at least four vertices. The edges and faces are ignored 
	 * @param numExtraPoints Number of extra random points added to the initial points
	 * @param glueEdges The edges 
	 * @param branchVertices
	 */
	public static void generateEllipticImage(CoHDS hds, int numExtraPoints, Set<CoEdge> glueEdges, int... branchVertices) {
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
		for (CoEdge e : new HashSet<CoEdge>(hds.getEdges())) {
			hds.removeEdge(e);
		}
		for (CoFace f : new HashSet<CoFace>(hds.getFaces())) {
			hds.removeFace(f);
		}
	
		// additional points
		Random rnd = new Random();
		for (int i = 0; i < numExtraPoints; i++) {
			double[] pos = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
			CoVertex v = hds.addNewVertex();
			v.setPosition(new Point(pos));
		}
		
		// on the sphere
		for (CoVertex v : hds.getVertices()) {
			Point vec = v.getPosition();
			vec.normalize();
		}
		
		// convex hull
		TypedAdapterSet<double[]> a = new TypedAdapterSet<double[]>(new CoPositionAdapter());
		ConvexHull.convexHull(hds, a, 1E-8);
		int vOffset = hds.numVertices();
		int eOffset = hds.numEdges();
		HalfEdgeUtils.copy(hds, hds);
		for (int i = 0; i < vOffset; i++) {
			CoVertex v = hds.getVertex(i);
			CoVertex vc = hds.getVertex(vOffset + i); 
			Point p = v.getPosition();
			Point p2 = new Point(p);
			vc.setPosition(p2);
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
		
		glueEdges.clear();
		glueEdges.addAll(path1);
		glueEdges.addAll(path2);
		glueEdges.addAll(path1c);
		glueEdges.addAll(path2c);
	}
	
	
	
	
	
	public static Complex calculateHalfPeriodRatioMathLink(Point p1, Point p2, Point p3, Point p4, KernelLink l) throws MathLinkException {
		// to the sphere
		p1.normalize();
		p2.normalize();
		p3.normalize();
		p4.normalize();
		// project stereographically
		Complex z1 = new Complex(p1.x() / (1 - p1.z()), p1.y() / (1 - p1.z()));
		Complex z2 = new Complex(p2.x() / (1 - p2.z()), p2.y() / (1 - p2.z()));
		Complex z3 = new Complex(p3.x() / (1 - p3.z()), p3.y() / (1 - p3.z()));
		Complex z4 = new Complex(p4.x() / (1 - p4.z()), p4.y() / (1 - p4.z()));
		Complex e1 = toInfinitZeroSum(z1, z1, z2, z3, z4);
		Complex e2 = toInfinitZeroSum(z2, z1, z2, z3, z4);
		Complex e3 = toInfinitZeroSum(z3, z1, z2, z3, z4);
		Complex g2 = e1.times(e2).plus(e2.times(e3)).plus(e3.times(e1)).times(-4);
		Complex g3 = e1.times(e2).times(e3).times(4);
		Complex[] w1w3 = invokeWeierstrassHalfperiods(g2, g3, l);
		return w1w3[1].divide(w1w3[0]);
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
		
		Point p1 = new Point(0.5257311122273196, 0.8506508082851774, -2.259622999145302E-9);
		Point p2 = new Point(-7.865783539640066E-9, -0.525731111594899, -0.8506508086760348);
		Point p3 = new Point(0.5257311164974198, -0.8506508056461103, -1.2759862959410077E-10);
		Point p4 = new Point(5.046344158026151E-10, 0.5257311116954958, 0.8506508086138624);
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
