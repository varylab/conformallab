package de.varylab.discreteconformal.unwrapper;

import geom3d.Point;

import java.util.Map;

import no.uib.cipr.matrix.Vector;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility;

public class SphericalUnwrapper implements Unwrapper{



	public Vector unwrap(CoHDS surface, int numCones, boolean quantizeCones) throws Exception {
		// punch out vertex 0 and reorder solver indices
		CoVertex v0 = surface.getVertex(0);
		for (CoEdge e : HalfEdgeUtils.incomingEdges(v0)) {
			Map<CoVertex, CoVertex> vMap = CuttingUtility.cutAtEdge(e);
			for (CoVertex vOld : vMap.keySet()) {
				CoVertex vNew = vMap.get(vOld);
				vNew.setPosition(vOld.getPosition());
			}
			if (numCones != 0) {
				break;
			}				
		}
		
//		diskUnwrapper.unwrap(hds);
		
//		HashSet<CVertex> boundary = new HashSet<CVertex>();
//		boundary.add(v0);
//		boundary.addAll(neighboringVertices(v0));
//		int n = hds.prepareInvariantData();
//		
//		// optimization
//		mon.subTask("Minimizing");
//		DenseVector u = new DenseVector(n);
//		Matrix H = new CompRowMatrix(n,n,makeNonZeros(hds));
//		NewtonOptimizer optimizer = new NewtonOptimizer(H);
//		optimizer.setStepController(new ArmijoStepController());
//		optimizer.setSolver(Solver.CG);
//		optimizer.setError(1E-5);
//		try {
//			optimizer.minimize(u, hds);
//		} catch (NotConvergentException e) {
//			mon.setCanceled(true);
//			throw new UnwrapException("Optimization did not succeed: " + e.getMessage());
//		}
//		mon.worked(1);
//		
//		// layout Euclidean
//		mon.subTask("Layout");
//		CDiskLayout.doLayout(hds, u, hds.calculateAlphas(u));
//		mon.worked(1);
//		
//		// spherical mapping
//		mon.subTask("Sphere mapping");
//		for (CVertex v : hds.getVertices()) {
//			Point t = v.getPosition();
//			t.setX(t.x() / t.z());
//			t.setY(t.y() / t.z());
//			t.setZ(1.0);
//		}
//		normalizeBeforeProjection(hds, 10.0);
//		inverseStereographicProjection(hds, 1.0);
//		v0.setTextureCoord(new Point(0.0, 0.0, 1.0));
//		try {
//			CSphereNormalizer.normalize(hds);
//		} catch (NotConvergentException e) {
//			mon.setCanceled(true);
//			throw new UnwrapException("Sphere normalization did not succeed: " + e.getMessage());
//		}
//		for (CVertex v : hds.getVertices()) {
//			v.setPosition(v.getTextureCoord());
//			Point p = v.getPosition();
//			Point t = v.getTextureCoord();
//			double U = Math.acos(p.z());
//			double V = Math.acos(p.x() / Math.sin(U));
//			t.setX(U);
//			t.setY(V);
//			t.setZ(1.0);
//		}
//		mon.worked(1);
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
	
	
	
	
}
