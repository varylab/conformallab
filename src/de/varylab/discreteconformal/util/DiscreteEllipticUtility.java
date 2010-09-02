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
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.calculator.SubdivisionCalculator;
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
			u = unwrapper.unwrap(hds);
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
		System.out.println("Using four vertices as branch points:");
		System.out.println(v1.getPosition());
		System.out.println(v2.getPosition());
		System.out.println(v3.getPosition());
		System.out.println(v4.getPosition());
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
		ConvexHull.convexHull(hds, new SubdivisionCalculator(), 1E-8);
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
		
		List<CoEdge> path1 = Search.bFS(v1, v2, new HashSet<CoVertex>());
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

	
}
