package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.unwrapper.BoundaryMode.Isometric;
import static de.varylab.discreteconformal.unwrapper.QuantizationMode.AllAngles;
import static de.varylab.discreteconformal.util.UnwrapUtility.prepareInvariantDataEuclidean;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import no.uib.cipr.matrix.DenseVector;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.jpetsc.InsertMode;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.ConvergenceFlags;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.MappedEdgeLengthAdapter;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanApplication;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanOptimizable;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;


public class StereographicUnwrapper implements Unwrapper{

	private Logger
		log = Logger.getLogger(StereographicUnwrapper.class.getName());
	private double
		gradTolerance = 1E-8;
	private int
		maxIterations = 150;
	private CoVertex
		unwrapRoot = null,
		layoutRoot = null;
	
	static {
		Tao.Initialize();
	}
	
	@Override
	public void unwrap(CoHDS surface, int genus, AdapterSet aSet) throws Exception {
		CoVertex v0 = unwrapRoot;
		if (v0 == null) {
			v0 = surface.getVertex(surface.numVertices() - 1);
		}
		// change edge lengths conformally such that the edges incident with vertex 0 are of equal length
		double meanLength = 0.0;
		int numIncident = 0;
		for (CoEdge e : HalfEdgeUtils.incomingEdges(v0)) {
			meanLength += aSet.get(Length.class, e, Double.class);
			numIncident++;
		}
		meanLength /= numIncident;
		Map<CoVertex, Double> uMap = new HashMap<CoVertex, Double>();
		for (CoEdge e : HalfEdgeUtils.incomingEdges(v0)) {
			double l = aSet.get(Length.class, e, Double.class);
			double u = meanLength / l;
			uMap.put(e.getStartVertex(), u);
		}
		Map<CoEdge, Double> lMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : HalfEdgeUtils.incomingEdges(v0)) {
			CoVertex v = e.getStartVertex();
			for (CoEdge ee : HalfEdgeUtils.incomingEdges(v)) {
				double l = aSet.get(Length.class, ee, Double.class);
				CoVertex vs = ee.getStartVertex();
				CoVertex vt = ee.getTargetVertex();
				double us = uMap.containsKey(vs) ? uMap.get(vs) : 1.0;
				double ut = uMap.containsKey(vt) ? uMap.get(vt) : 1.0;
				double nl = l * us * ut;
				lMap.put(ee, nl);
				lMap.put(ee.getOppositeEdge(), nl);
				System.out.println(ee + ": " + l + " -> " + nl);
			}
		}
		MappedEdgeLengthAdapter eAdapter = new MappedEdgeLengthAdapter(lMap, 100);
		aSet.add(eAdapter);
		
		// punch out the last vertex
		TopologyAlgorithms.removeVertex(v0);
		CEuclideanOptimizable opt = new CEuclideanOptimizable(surface);
		prepareInvariantDataEuclidean(opt.getFunctional(), surface, Isometric, AllAngles, aSet);
		
		// optimization
//		DenseVector u = new DenseVector(n);
//		Matrix H = opt.getHessianTemplate();
//		NewtonOptimizer optimizer = new NewtonOptimizer(H);
//		optimizer.setStepController(new ArmijoStepController());
//		optimizer.setSolver(Solver.BiCGstab);
//		optimizer.setError(gradTolerance);
//		optimizer.setMaxIterations(maxIterations);
//		try {
//			optimizer.minimize(u, opt);
//		} catch (NotConvergentException e) {
//			throw new UnwrapException("Optimization did not succeed: " + e.getMessage());
//		}
		
		CEuclideanApplication app = new CEuclideanApplication(surface);
		int n = app.getDomainDimension();
		Vec u = new Vec(n);
		// set variable lambda start values
		for (CoEdge e : surface.getPositiveEdges()) {
			if (e.getSolverIndex() >= 0) {
				u.setValue(e.getSolverIndex(), e.getLambda(), InsertMode.INSERT_VALUES);
			}
		}
		app.setInitialSolutionVec(u);
		Mat H = app.getHessianTemplate();
		app.setHessianMat(H, H);
		
		Tao optimizer = new Tao(Tao.Method.NTR);
		optimizer.setApplication(app);
		optimizer.setGradientTolerances(gradTolerance, gradTolerance, gradTolerance); 
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setMaximumIterates(maxIterations);
		System.out.println("Using grad tolerance " + gradTolerance);
		optimizer.solve();
		if (optimizer.getSolutionStatus().reason != ConvergenceFlags.CONVERGED_ATOL) {
			throw new RuntimeException("Optinizer did not converge: \n" + optimizer.getSolutionStatus());
		}
		System.out.println(optimizer.getSolutionStatus());
		double[] uValues = u.getArray();
		DenseVector uVec = new DenseVector(uValues);
		
		// layout Euclidean
		layoutRoot = EuclideanLayout.doLayout(surface, opt.getFunctional(), uVec);
//		if (true) return;
		// spherical mapping
		for (CoVertex v : surface.getVertices()) {
			Pn.dehomogenize(v.T, v.T);
		}
		normalizeBeforeProjection(surface, 1);
		inverseStereographicProjection(surface, 1.0);
		
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
		assert firstPrev != null;
		firstPrev.linkOppositeEdge(lastNext);
		
		try {
			SphericalNormalizerPETSc.normalize(surface, aSet, TexturePosition4d.class, TexturePosition.class);
		} catch (Exception e) {
			log.info("Sphere normalization did not succeed: " + e.getMessage());
		}
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
		unwrapRoot = root;
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
