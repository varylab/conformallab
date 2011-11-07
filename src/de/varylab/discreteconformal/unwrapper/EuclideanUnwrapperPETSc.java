package de.varylab.discreteconformal.unwrapper;

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;

import no.uib.cipr.matrix.DenseVector;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.jpetsc.InsertMode;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.jtem.jtao.Tao.GetSolutionStatusResult;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanApplication;
import de.varylab.discreteconformal.unwrapper.numerics.MTJDomain;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.Search.DefaultWeightAdapter;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.UnwrapUtility.BoundaryMode;
import de.varylab.discreteconformal.util.UnwrapUtility.QuantizationMode;

public class EuclideanUnwrapperPETSc implements Unwrapper {

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

	public static double
		lastGNorm = 0;
	private Collection<CoVertex>
		cones = new HashSet<CoVertex>();
	
	public CoVertex
		layoutRoot = null;
	public CuttingInfo<CoVertex, CoEdge, CoFace> 
		cutInfo = null;
	public Map<CoEdge, Double>
		lengthMap = null;
	
	static {
		Tao.Initialize();		
	}
	
	@Override
	public void unwrap(CoHDS surface, int genus, AdapterSet aSet) throws Exception {
		CEuclideanApplication app = new CEuclideanApplication(surface);
		double[] uValues = calculateConformalFactors(surface, aSet, app);
		DenseVector uVec = new DenseVector(uValues);

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
		lengthMap = EuclideanLayout.getLengthMap(surface, app.getFunctional(), uVec);
		layoutRoot = EuclideanLayout.doLayout(surface, app.getFunctional(), uVec);
	}

	
	public Map<CoVertex, Double> calculateConformalFactors(CoHDS surface, AdapterSet aSet) throws UnwrapException {
		CEuclideanApplication app = new CEuclideanApplication(surface);
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
	
	
	private double[] calculateConformalFactors(CoHDS surface, AdapterSet aSet, CEuclideanApplication app) throws UnwrapException {
		UnwrapUtility.prepareInvariantDataEuclidean(app.getFunctional(), surface, boundaryMode, boundaryQuantMode, aSet);
		// cones
		cones = ConesUtility.setUpCones(surface, numCones); 
		// optimization
		Vec u;
		Tao optimizer;
		int n = app.getDomainDimension();
		u = new Vec(n);
		// set variable lambda start values
		for (CoEdge e : surface.getPositiveEdges()) {
			if (e.getSolverIndex() >= 0) {
				u.setValue(e.getSolverIndex(), e.getLambda(), InsertMode.INSERT_VALUES);
			}
		}
		app.setInitialSolutionVec(u);
		Mat H = app.getHessianTemplate();
		app.setHessianMat(H, H);

		optimizer = new Tao(Tao.Method.NTR);
		optimizer.setApplication(app);
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setGradientTolerances(gradTolerance, gradTolerance, gradTolerance);
		optimizer.setMaximumIterates(maxIterations);
		optimizer.solve();
		
		GetSolutionStatusResult status = optimizer.getSolutionStatus();
		lastGNorm = status.gnorm;
		if (status.reason.cvalue() < 0) {
			throw new UnwrapException("Optimization did not succeed: " + status);
		}

		if (!cones.isEmpty()) {
			if (conesMode != QuantizationMode.AllAngles) {
				// calculating cones
				cones = ConesUtility.quantizeCones(surface, cones, conesMode);
				
				// optimizing conformal structure
				CEuclideanApplication app2 = new CEuclideanApplication(surface);
				n = app2.getDomainDimension();
	
				u = new Vec(n);
				H = app.getHessianTemplate();
				
				app2.setInitialSolutionVec(u);
				app2.setHessianMat(H, H);
				
				optimizer.setApplication(app2);
				optimizer.setGradientTolerances(1E-5, 0, 0);
				optimizer.solve();
				
				status = optimizer.getSolutionStatus();
				System.out.println("Cone Quantization: " + status);
				if (status.reason.cvalue() < 0) {
					throw new UnwrapException("Cone quantization did not succeed: " + status);
				}
			}
		}
		double [] uValues = u.getArray();
		u.restoreArray();
		return uValues;
	}
	
	public Collection<CoVertex> getCones() {
		return cones;
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
