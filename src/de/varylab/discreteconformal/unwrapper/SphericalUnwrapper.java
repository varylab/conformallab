package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.util.UnwrapUtility.prepareInvariantDataEuclidean;
import static de.varylab.discreteconformal.util.UnwrapUtility.BoundaryMode.Isometric;
import static de.varylab.discreteconformal.util.UnwrapUtility.QuantizationMode.AllAngles;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanOptimizable;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;
import de.varylab.mtjoptimization.stepcontrol.ArmijoStepController;


public class SphericalUnwrapper implements Unwrapper{

	private Logger
		log = Logger.getLogger(SphericalUnwrapper.class.getName());
	private double
		gradTolerance = 1E-8;
	private int
		maxIterations = 150;
	private CoVertex
		layoutRoot = null;

	
	@Override
	public void unwrap(CoHDS surface, int genus, AdapterSet aSet) throws Exception {
		// change edge lengths conformally such that the edges incident with vertex 0 are of equal length
		
		
		
		// punch out the last vertex
		CoVertex v0 = surface.getVertex(surface.numVertices() - 1);
		TopologyAlgorithms.removeVertex(v0);
		CEuclideanOptimizable opt = new CEuclideanOptimizable(surface);
		int n = prepareInvariantDataEuclidean(opt.getFunctional(), surface, Isometric, AllAngles, aSet);
		
		// optimization
		DenseVector u = new DenseVector(n);
		Matrix H = opt.getHessianTemplate();
		NewtonOptimizer optimizer = new NewtonOptimizer(H);
		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.BiCGstab);
		optimizer.setError(gradTolerance);
		optimizer.setMaxIterations(maxIterations);
		
		try {
			optimizer.minimize(u, opt);
		} catch (NotConvergentException e) {
			throw new UnwrapException("Optimization did not succeed: " + e.getMessage());
		}
		
		// layout Euclidean
		layoutRoot = EuclideanLayout.doLayout(surface, opt.getFunctional(), u);
		
		// spherical mapping
		for (CoVertex v : surface.getVertices()) {
			Pn.dehomogenize(v.T, v.T);
		}
		normalizeBeforeProjection(surface, 1);
		inverseStereographicProjection(surface, 1.0);
		try {
			SphericalNormalizer.normalize(surface);
		} catch (NotConvergentException e) {
			log.info("Sphere normalization did not succeed: " + e.getMessage());
		}
		
		// re-insert vertex 0
		CoVertex oldV0 = v0;
		v0 = surface.addNewVertex();
		v0.P = oldV0.P;
		v0.T = new double[] {0,1,0,1};
		
		CoEdge lastNext = null;
		CoEdge firstPrev = null;
		Set<CoEdge> boundary = new HashSet<CoEdge>(HalfEdgeUtils.boundaryEdges(surface));
		CoEdge be = HalfEdgeUtils.boundaryEdges(surface).iterator().next();
		while (!boundary.isEmpty()) {
			CoEdge next = be.getNextEdge();
			boundary.remove(be); 
			CoEdge eNext = surface.addNewEdge();
			CoEdge ePrev = surface.addNewEdge();
			CoFace f = surface.addNewFace();
			be.linkNextEdge(eNext);
			be.linkPreviousEdge(ePrev);
			eNext.linkNextEdge(ePrev);
			eNext.setTargetVertex(v0);
			ePrev.setTargetVertex(be.getStartVertex());
			be.setLeftFace(f);
			eNext.setLeftFace(f);
			ePrev.setLeftFace(f);
			if (lastNext != null) {
				ePrev.linkOppositeEdge(lastNext);
			}
			if (firstPrev == null) {
				firstPrev = ePrev;
			}
			lastNext = eNext;
			be = next;
		}
		firstPrev.linkOppositeEdge(lastNext);
	}

	
	/**
	 * Project stereographically onto the sphere
	 * @param graph
	 * @param scale
	 */
	public static void inverseStereographicProjection(CoHDS hds, double scale){
		for (CoVertex v : hds.getVertices()){
			Pn.dehomogenize(v.T, v.T);
			double x = v.T[0] / scale;
			double y = v.T[1] / scale;
			double nx = 2 * x;
			double ny = x*x + y*y - 1;
			double nz = 2 * y;
			double nw = ny + 2;
			v.T = new double[] {nx / nw, ny / nw, nz / nw, 1.0};
		}
	}
	
	
	public static double[] baryCenter(CoHDS hds){
		double[] result = {0,0,0,1};
		for (CoVertex v : hds.getVertices()){
			Rn.add(result, v.T, result);
		}
		return Rn.times(result, 1.0 / hds.numVertices(), result);
	}
	
	
	public static double meanRadius(CoHDS hds){
		double result = 0;
		double[] zero = {0,0,0,1};
		for (CoVertex v : hds.getVertices()){
			result += Pn.distanceBetween(zero, v.T, Pn.EUCLIDEAN);
		}
		return result / hds.numVertices();
	}
	
	public static void normalizeBeforeProjection(CoHDS hds, double scale){
		double[] offset = baryCenter(hds);
		Pn.dehomogenize(offset, offset);
		for (CoVertex v : hds.getVertices()){
			Rn.subtract(v.T, v.T, offset);
			v.T[3] = 1.0;
		}
		scale = meanRadius(hds) / scale;
		for (CoVertex v : hds.getVertices()){
			v.T[3] *= scale;
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

	@Override
	public void setCutRoot(CoVertex root) {
		
	}

	@Override
	public CuttingInfo<CoVertex, CoEdge, CoFace> getCutInfo() {
		return new CuttingInfo<CoVertex, CoEdge, CoFace>();
	}
	@Override
	public Map<CoEdge, Double> getlengthMap() {
		return null;
	}
	@Override
	public CoVertex getLayoutRoot() {
		return layoutRoot;
	}
	
}
