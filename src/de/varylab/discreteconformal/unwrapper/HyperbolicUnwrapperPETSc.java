package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.util.CuttingUtility.cutManifoldToDisk;
import static de.varylab.discreteconformal.util.SparseUtility.getPETScNonZeros;

import java.util.Map;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.jpetsc.InsertMode;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.PETSc;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.ConvergenceFlags;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.adapter.HyperbolicLengthWeightAdapter;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperbolicApplication;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.UnwrapUtility;

public class HyperbolicUnwrapperPETSc implements Unwrapper {

	private double
		gradTolerance = 1E-8;
	private int
		maxIterations = 150;
	private CoVertex
		cutRoot = null;
	private boolean
		cutAndLayout = true;
	
	private CoVertex
		layoutRoot = null;
	private CuttingInfo<CoVertex, CoEdge, CoFace> 
		cutInfo = null;
	private Map<CoEdge, Double>
		lengthMap = null;
	private Vector
		uVec = null;
	
	static {
		Tao.Initialize();		
	}
	
	@Override
	public void unwrap(CoHDS surface, int genus, AdapterSet aSet) throws Exception {
		CHyperbolicApplication app = new CHyperbolicApplication(surface);
		int n = app.getDomainDimension(); 
		Vec u = new Vec(n);
		
		UnwrapUtility.prepareInvariantDataHyperbolic(app.getFunctional(), surface, aSet);
		
		// set variable lambda start values
		boolean hasCircularEdges = false;
		for (CoEdge e : surface.getPositiveEdges()) {
			if (e.getSolverIndex() >= 0) {
				u.setValue(e.getSolverIndex(), e.getLambda(), InsertMode.INSERT_VALUES);
				hasCircularEdges = true;
			}
		}
		app.setInitialSolutionVec(u);
		if (!hasCircularEdges) {
			Mat H = Mat.createSeqAIJ(n, n, PETSc.PETSC_DEFAULT, getPETScNonZeros(surface));
			H.assemble();
			app.setHessianMat(H, H);
		}
		
		Tao optimizer = new Tao(hasCircularEdges ? Tao.Method.LMVM : Tao.Method.NTR);
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
		uVec = new DenseVector(u.getArray());
		u.restoreArray();

		if (cutAndLayout) {
			HyperbolicLengthWeightAdapter hypWa = new HyperbolicLengthWeightAdapter(uVec);
			if (cutRoot == null) {
				cutRoot = surface.getVertex(HyperbolicUnwrapper.getMinUIndex(uVec));
			}
			cutInfo = cutManifoldToDisk(surface, cutRoot, hypWa);
			lengthMap = HyperbolicLayout.getLengthMap(surface, uVec);
			layoutRoot = HyperbolicLayout.doLayout(surface, null, uVec);
		}
	}

	
	public Vector getUResult() {
		return uVec;
	}
	
	public boolean isCutAndLayout() {
		return cutAndLayout;
	}
	public void setCutAndLayout(boolean cutAndLayout) {
		this.cutAndLayout = cutAndLayout;
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
	
}
