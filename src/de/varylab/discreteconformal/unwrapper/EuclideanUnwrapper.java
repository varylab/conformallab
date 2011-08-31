package de.varylab.discreteconformal.unwrapper;

import java.util.Collection;
import java.util.Map;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanOptimizable;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.Search.DefaultWeightAdapter;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.UnwrapUtility.BoundaryMode;
import de.varylab.discreteconformal.util.UnwrapUtility.QuantizationMode;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;
import de.varylab.mtjoptimization.stepcontrol.ArmijoStepController;

public class EuclideanUnwrapper implements Unwrapper {

	private QuantizationMode
		conesMode = QuantizationMode.AllAngles,
		boundaryQuantMode = QuantizationMode.AllAngles;
	private BoundaryMode
		boundaryMode = BoundaryMode.Isometric;
	private int 
		maxIterations = 150,
		numCones = 0;
	private double
		gradTolerance = 1E-8;
	
	public CoVertex
		layoutRoot = null;
	public CuttingInfo<CoVertex, CoEdge, CoFace> 
		cutInfo = null;
	public Map<CoEdge, Double>
		lengthMap = null;
	
	@Override
	public void unwrap(CoHDS surface, int genus, AdapterSet aSet) throws Exception {
		CEuclideanOptimizable opt = new CEuclideanOptimizable(surface);
		UnwrapUtility.prepareInvariantDataEuclidean(opt.getFunctional(), surface, boundaryMode, boundaryQuantMode, aSet);
		// cones
		Collection<CoVertex> cones = ConesUtility.setUpCones(surface, numCones); 
		// optimization
		int n = opt.getDomainDimension();
		DenseVector u = new DenseVector(n);
		// set variable lambda start values
		for (CoEdge e : surface.getPositiveEdges()) {
			if (e.getSolverIndex() >= 0) {
				u.set(e.getSolverIndex(), e.getLambda());
			}
		}
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
		if (!cones.isEmpty()) {
			if (conesMode != QuantizationMode.AllAngles) {
				cones = ConesUtility.quantizeCones(surface, cones, conesMode);
				n = opt.getDomainDimension();
				u = new DenseVector(n);
				H = opt.getHessianTemplate();
				optimizer.setHessianTemplate(H);
				try {
					optimizer.minimize(u, opt);
				} catch (NotConvergentException e) {
					throw new UnwrapException("Cone quantization did not succeed: " + e.getMessage());
				}
			}
		}
		
		switch (genus) {
		case 0:
			cutInfo = ConesUtility.cutMesh(surface);
			break;
		case 1:
			CoVertex cutRoot = surface.getVertex(0);
			DefaultWeightAdapter<CoEdge> constantWeight = new DefaultWeightAdapter<CoEdge>();
			cutInfo = CuttingUtility.cutTorusToDisk(surface, cutRoot, constantWeight);
			break;
		default:
			throw new RuntimeException("Cannot work on higher genus");
		}
		lengthMap = EuclideanLayout.getLengthMap(surface, opt.getFunctional(), u);
		layoutRoot = EuclideanLayout.doLayout(surface, opt.getFunctional(), u);
	}
	
	public void setNumCones(int numCones) {
		this.numCones = numCones;
	}
	public void setConeMode(QuantizationMode quantizationMode) {
		this.conesMode = quantizationMode;
	}
	@Override
	public void setGradientTolerance(double tol) {
		gradTolerance = tol;
	}
	@Override
	public void setMaxIterations(int maxIterations) {
		this.maxIterations = maxIterations;
	}
	public void setBoundaryQuantMode(QuantizationMode boundaryQuantMode) {
		this.boundaryQuantMode = boundaryQuantMode;
	}
	public void setBoundaryMode(BoundaryMode boundaryMode) {
		this.boundaryMode = boundaryMode;
	}
	@Override
	public void setCutRoot(CoVertex root) {
	}

	@Override
	public CuttingInfo<CoVertex, CoEdge, CoFace> getCutInfo() {
		return cutInfo;
	}
	@Override
	public Map<CoEdge, Double> getlengthMap() {
		return lengthMap;
	}
	@Override
	public CoVertex getLayoutRoot() {
		return layoutRoot;
	}
	
}
