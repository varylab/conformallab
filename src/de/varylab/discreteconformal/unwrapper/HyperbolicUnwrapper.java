package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.util.CuttingUtility.cutManifoldToDisk;

import java.util.Map;
import java.util.Set;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.varylab.discreteconformal.adapter.HyperbolicLengthWeightAdapter;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperbolicOptimizable;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.UnwrapUtility.ZeroU;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;
import de.varylab.mtjoptimization.stepcontrol.ArmijoStepController;

public class HyperbolicUnwrapper implements Unwrapper {

	private double
		gradTolerance = 1E-8;
	private int
		maxIterations = 150;
	private CoVertex
		cutRoot = null;
	
	public CoVertex
		layoutRoot = null;
	public CuttingInfo<CoVertex, CoEdge, CoFace> 
		cutInfo = null;
	public Map<CoEdge, Double>
		lengthMap = null;
	public Set<CoEdge>
		cutGraph = null;
	
	
	@Override
	public void unwrap(CoHDS surface, int genus, AdapterSet aSet) throws Exception {
		CHyperbolicOptimizable opt = new CHyperbolicOptimizable(surface);
		ZeroU zeroU = new ZeroU();
		UnwrapUtility.prepareInvariantDataHyperbolicAndSpherical(opt.getFunctional(), surface, aSet, zeroU);
		
		// optimization
		int n = opt.getDomainDimension();
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

		HyperbolicLengthWeightAdapter hypWa = new HyperbolicLengthWeightAdapter(u);
		if (cutRoot == null) {
			cutRoot = surface.getVertex(getMinVertexUIndex(u, surface.numVertices()));
		}
		if (cutGraph != null) {
			cutInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
			cutInfo.cutRoot = cutRoot;
			CuttingUtility.cutAtEdges(cutInfo, cutGraph);
		} else {
			cutInfo = cutManifoldToDisk(surface, cutRoot, hypWa);
		}
		lengthMap = HyperbolicLayout.getLengthMap(surface, opt.getFunctional(), u);
		layoutRoot = HyperbolicLayout.doLayout(surface, null, opt.getFunctional(), u);
	}
	

	protected static int getMinVertexUIndex(Vector u, int numVertices) {
		int index = 0;
		double iVal = u.get(0);
		for (int i = 1; i < Math.min(u.size(), numVertices); i++) {
			double val = u.get(i);
			if (iVal < val) {
				index = i;
				iVal = i;
			}
		}
		return index;
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
		this.cutRoot = root;
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
	@Override
	public void setCutGraph(Set<CoEdge> cutEdges) {
		this.cutGraph = cutEdges;
	}

}
