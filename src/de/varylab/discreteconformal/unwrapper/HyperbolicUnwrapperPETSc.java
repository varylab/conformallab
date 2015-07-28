package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.unwrapper.HyperbolicUnwrapper.getMinVertexUIndex;
import static de.varylab.discreteconformal.util.CuttingUtility.cutManifoldToDisk;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.ConvergenceFlags;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.adapter.HyperbolicLengthWeightAdapter;
import de.varylab.discreteconformal.functional.ConformalFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperbolicApplication;
import de.varylab.discreteconformal.unwrapper.numerics.MTJDomain;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.UnwrapUtility.ZeroU;

public class HyperbolicUnwrapperPETSc implements Unwrapper {

	private Logger
		log = Logger.getLogger(getClass().getName());
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
	private Set<CoEdge>
		cutGraph = null;
	private CuttingInfo<CoVertex, CoEdge, CoFace> 
		cutInfo = null;
	private Map<CoEdge, Double>
		lengthMap = null;
	private Vector
		uVec = null;
	private CHyperbolicApplication
		app = null;
	
	static {
		Tao.Initialize();		
	}
	
	public HyperbolicUnwrapperPETSc() {
		super();
	}

	public HyperbolicUnwrapperPETSc(boolean cutAndLayout) {
		super();
		this.cutAndLayout = cutAndLayout;
	}

	@Override
	public void unwrap(CoHDS surface, int genus, AdapterSet aSet) throws Exception {
		app = new CHyperbolicApplication(surface);
		double[] uArr = calculateConformalFactors(surface, aSet, app);
		uVec = new DenseVector(uArr);

		if (cutAndLayout) {
			HyperbolicLengthWeightAdapter hypWa = new HyperbolicLengthWeightAdapter(uVec);
			if (cutRoot == null) {
				cutRoot = surface.getVertex(getMinVertexUIndex(uVec, surface.numVertices()));
			}
			if (cutGraph != null) {
				cutInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
				cutInfo.cutRoot = cutRoot;
				CuttingUtility.cutAtEdges(cutInfo, cutGraph);
			} else {
				cutInfo = cutManifoldToDisk(surface, cutRoot, hypWa);
			}
			lengthMap = HyperbolicLayout.getLengthMap(surface, app.getFunctional(), uVec);
			layoutRoot = HyperbolicLayout.doLayout(surface, null, app.getFunctional(), uVec);
		}
	}
	
	
	public Map<CoVertex, Double> calculateConformalFactors(CoHDS surface, AdapterSet aSet) throws UnwrapException {
		app = new CHyperbolicApplication(surface);
		double[] u = calculateConformalFactors(surface, aSet, app);
		DenseVector uVec = new DenseVector(u);
		MTJDomain uDomain = new MTJDomain(uVec);
		Map<CoVertex, Double> result = new LinkedHashMap<CoVertex, Double>();
		for (CoVertex v : surface.getVertices()) {
			Double uVal = app.getFunctional().getVertexU(v, uDomain);
			result.put(v, uVal);
		}
		return result;
	}


	private double[] calculateConformalFactors(CoHDS surface, AdapterSet aSet, CHyperbolicApplication app) {
		ZeroU zeroU = new ZeroU();
		UnwrapUtility.prepareInvariantDataHyperbolicAndSpherical(app.getFunctional(), surface, aSet, zeroU);
		int n = app.getDomainDimension(); 
		Vec u = new Vec(n);
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
		UnwrapUtility.logSolutionStatus(optimizer, log);
		double[] uVec = u.getArray();
		u.restoreArray();
		return uVec;
	}

	
	public ConformalFunctional<CoVertex, CoEdge, CoFace> getFunctional() {
		return app.getFunctional();
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
	@Override
	public void setCutGraph(Set<CoEdge> cutEdges) {
		this.cutGraph = cutEdges;
	}
	
}
