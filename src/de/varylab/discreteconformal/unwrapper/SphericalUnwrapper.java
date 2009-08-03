package de.varylab.discreteconformal.unwrapper;

import geom3d.Point;

import java.util.Map;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanOptimizable;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.SparseUtility;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;
import de.varylab.mtjoptimization.stepcontrol.ArmijoStepController;

public class SphericalUnwrapper implements Unwrapper{

	private double
		gradTolerance = 1E-8;
	private int
		maxIterations = 150;

	/**
	 * This is completely untested and probably rubbish
	 */
	public Vector unwrap(CoHDS surface) throws Exception {
		// punch out vertex 0 and reorder solver indices
		CoVertex v0 = surface.getVertex(0);
		for (CoEdge e : HalfEdgeUtils.incomingEdges(v0)) {
			Map<CoVertex, CoVertex> vMap = CuttingUtility.cutAtEdge(e);
			for (CoVertex vOld : vMap.keySet()) {
				CoVertex vNew = vMap.get(vOld);
				vNew.setPosition(vOld.getPosition());
			}			
		}
		
//		diskUnwrapper.unwrap(hds);
		
//		HashSet<CoVertex> boundary = new HashSet<CoVertex>();
//		boundary.add(v0);
//		boundary.addAll(neighboringVertices(v0));
		int n = surface.prepareInvariantDataEuclidean();
		CEuclideanOptimizable opt = new CEuclideanOptimizable(surface);
		
		// optimization
		DenseVector u = new DenseVector(n);
		Matrix H = new CompRowMatrix(n,n, SparseUtility.makeNonZeros(surface));
		NewtonOptimizer optimizer = new NewtonOptimizer(H);
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.CG);
		optimizer.setError(gradTolerance);
		optimizer.setMaxIterations(maxIterations);
		try {
			optimizer.minimize(u, opt);
		} catch (NotConvergentException e) {
			throw new UnwrapException("Optimization did not succeed: " + e.getMessage());
		}
		
		// layout Euclidean
		EuclideanLayout.doLayout(surface, u);
		
		// spherical mapping
		for (CoVertex v : surface.getVertices()) {
			Point t = v.getPosition();
			t.setX(t.x() / t.z());
			t.setY(t.y() / t.z());
			t.setZ(1.0);
		}
		normalizeBeforeProjection(surface, 10.0);
		inverseStereographicProjection(surface, 1.0);
		v0.setTextureCoord(new Point(0.0, 0.0, 1.0));
		try {
			SphericalNormalizer.normalize(surface);
		} catch (NotConvergentException e) {
			throw new UnwrapException("Sphere normalization did not succeed: " + e.getMessage());
		}
		for (CoVertex v : surface.getVertices()) {
			v.setPosition(v.getTextureCoord());
			Point p = v.getPosition();
			Point t = v.getTextureCoord();
			double U = Math.acos(p.z());
			double V = Math.acos(p.x() / Math.sin(U));
			t.setX(U);
			t.setY(V);
			t.setZ(1.0);
		}
		return null;
	}

	
	
	/**
	 * Project stereographically onto the sphere
	 * @param graph
	 * @param scale
	 */
	public static void inverseStereographicProjection(CoHDS hds, double scale){
		for (CoVertex v : hds.getVertices()){
			double x = v.getTextureCoord().x() / scale;
			double y = v.getTextureCoord().y() / scale;
			double nx = 2 * x;
			double ny = x*x + y*y - 1;
			double nz = 2 * y;
			double nw = ny + 2;
			v.getTextureCoord().set(nx / nw, ny / nw, nz / nw);
		}
	}
	
	
	public static Point baryCenter(CoHDS hds){
		Point result = new Point(0,0,0);
		for (CoVertex v : hds.getVertices()){
			result.add(v.getTextureCoord());
		}
		result.times(1.0 / hds.numVertices());
		return result;
	}
	
	
	public static double meanRadius(CoHDS hds){
		double result = 0;
		for (CoVertex v : hds.getVertices()){
			result += v.getTextureCoord().getLength();
		}
		return result / hds.numVertices();
	}
	
	public static void normalizeBeforeProjection(CoHDS hds, double scale){
		Point offset = baryCenter(hds);
		for (CoVertex v : hds.getVertices()){
			v.getTextureCoord().subtract(offset);
		}
		scale = meanRadius(hds) / scale;
		for (CoVertex v : hds.getVertices()){
			v.getTextureCoord().times(1 / scale);
		}
	
	}
	

	
	@Override
	public void setGradientTolerance(double tol) {
		gradTolerance = tol;
	}
	
	@Override
	public void setMaxIterations(int maxIterations) {
		this.maxIterations = maxIterations;
	}
	
	
	
}
