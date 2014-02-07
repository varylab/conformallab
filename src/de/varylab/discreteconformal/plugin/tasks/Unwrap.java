package de.varylab.discreteconformal.plugin.tasks;

import static java.lang.Math.PI;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import de.jreality.plugin.job.AbstractJob;
import de.jtem.blas.ComplexMatrix;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.theta.SiegelReduction;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.BoundaryMode;
import de.varylab.discreteconformal.unwrapper.CircleDomainUnwrapper;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapper;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.HyperbolicUnwrapper;
import de.varylab.discreteconformal.unwrapper.HyperbolicUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.QuantizationMode;
import de.varylab.discreteconformal.unwrapper.StereographicUnwrapper;
import de.varylab.discreteconformal.unwrapper.Unwrapper;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.DiscreteEllipticUtility;

public class Unwrap extends AbstractJob {

	private CoHDS
		surface = null;
	private AdapterSet
		aSet = new AdapterSet();
	private Set<CoVertex>
		selectedVertices = new HashSet<CoVertex>();
	private Set<CoEdge>
		selectedEdges = new HashSet<CoEdge>();
	private boolean
		usePetsc = false,
		spherize = false,
		useSelectionCuts = false;
	private QuantizationMode
		coneMode = QuantizationMode.AllAngles,
		boundaryQuantMode = QuantizationMode.AllAngles;
	private BoundaryMode
		boundaryMode = BoundaryMode.Isometric;
	private int
		numCones = 0,
		maxIterations = 150;
	private double
		toleranceExp = -8;
	
	// result values -------
	public int 
		genus = -1;
	public CoVertex
		layoutRoot = null;
	public CuttingInfo<CoVertex, CoEdge, CoFace> 
		cutInfo = null;
	public Map<CoEdge, Double>
		lengthMap = null;
	
	
	
	public Unwrap(CoHDS surface, AdapterSet aSet) {
		this.surface = surface;
		this.aSet = aSet;
	}
	
	@Override
	public String getJobName() {
		return "Discrete Conformal";
	}
	
	
	@Override
	public void executeJob() throws Exception {
		long startTime = System.currentTimeMillis();
		long unwrapTime = -1;
		long layoutTime = -1;
		fireJobProgress(0.0);
		if (!HalfEdgeUtils.isValidSurface(surface, true)) {
			throw new RuntimeException("Surface is not valid");
		}
		genus = HalfEdgeUtils.getGenus(surface);
		if (spherize) genus = 0;
		double gradTolerance = Math.pow(10, toleranceExp);
		switch (genus) {
		// disk or sphere---------------------
		case 0: 
			Collection<CoEdge> bList = HalfEdgeUtils.boundaryEdges(surface);
			boolean isSphere = bList.size() == 0;
			Unwrapper unwrapper = null;
			if (isSphere) {
				System.out.println("unwrapping a sphere...");
				if (usePetsc) {
//					unwrapper = new SphericalUnwrapperPETSc();
					unwrapper = new StereographicUnwrapper();
				} else {
//					unwrapper = new SphericalUnwrapper();
					unwrapper = new StereographicUnwrapper();
				}
			} else {
				if (boundaryMode == BoundaryMode.Circle) {
					CircleDomainUnwrapper.flipAtEars(surface, aSet);
					CircleDomainUnwrapper cdu = new CircleDomainUnwrapper();
					unwrapper = cdu;
				} else {
					if (usePetsc) {
						EuclideanUnwrapperPETSc uw = new EuclideanUnwrapperPETSc();
						uw.setNumCones(numCones);
						uw.setConeMode(coneMode);
						uw.setBoundaryMode(boundaryMode);
						uw.setBoundaryQuantMode(boundaryQuantMode);
						unwrapper = uw;
					} else {
						EuclideanUnwrapper uw = new EuclideanUnwrapper();
						uw.setNumCones(numCones);
						uw.setConeMode(coneMode);
						uw.setBoundaryMode(boundaryMode);
						uw.setBoundaryQuantMode(boundaryQuantMode);
						unwrapper = uw;
					}
				}
			}
			if (!selectedVertices.isEmpty()) {
				CoVertex cutRoot = selectedVertices.iterator().next();
				unwrapper.setCutRoot(cutRoot);
			}
			unwrapper.setGradientTolerance(gradTolerance);
			unwrapper.setMaxIterations(maxIterations);
			unwrapper.unwrap(surface, genus, aSet);
			unwrapTime = System.currentTimeMillis();
			fireJobProgress(0.5);
			
			// check sphere conditions
			for (CoFace f : surface.getFaces()) {
				CoEdge a = f.getBoundaryEdge();
				CoEdge b = a.getNextEdge();
				CoEdge c = b.getNextEdge();
				double A = a.getAlpha() + b.getAlpha() + c.getAlpha() - PI;
				if (A > 2*PI) {
					System.err.println("face " + f + " has an area greater than 2PI");
				}
			}
//			for (CoVertex v : surface.getVertices()) {
//				double sum = 0.0;
//				for (CoEdge e : HalfEdgeUtils.incomingEdges(v)) {
//					sum += e.getPreviousEdge().getAlpha();
//				}
//				if (!HalfEdgeUtils.isBoundaryVertex(v) && Math.abs(sum - v.getTheta()) > 1E-5) {
//					System.err.println("angle sum at vertex " + v + " is incorrect: expected " + v.getTheta() + ", actual: " + sum);
//				}
//			}
			
			cutInfo = unwrapper.getCutInfo();
			lengthMap = unwrapper.getlengthMap();
			layoutRoot = unwrapper.getLayoutRoot();
			layoutTime = System.currentTimeMillis();
			fireJobProgress(1.0);
			break;
		// torus ----------------------------
		case 1:
			System.out.println("unwrapping a torus...");
			if (usePetsc) {
				unwrapper = new EuclideanUnwrapperPETSc();
			} else {
				unwrapper = new EuclideanUnwrapper();
			}
			if (useSelectionCuts) unwrapper.setCutGraph(selectedEdges);
			if (!selectedVertices.isEmpty()) {
				CoVertex cutRoot = selectedVertices.iterator().next();
				unwrapper.setCutRoot(cutRoot);
			}
			unwrapper.setGradientTolerance(gradTolerance);
			unwrapper.setMaxIterations(maxIterations);
			unwrapper.unwrap(surface, genus, aSet);
			unwrapTime = System.currentTimeMillis();
			fireJobProgress(0.5);
			
			cutInfo = unwrapper.getCutInfo();
			lengthMap = unwrapper.getlengthMap();
			layoutRoot = unwrapper.getLayoutRoot();
			layoutTime = System.currentTimeMillis();
			try {
				Complex tau = DiscreteEllipticUtility.calculateHalfPeriodRatio(cutInfo);
				Complex piNormalizedTau = new Complex(0, 2 * PI);
				piNormalizedTau = piNormalizedTau.times(tau);
				final ComplexMatrix PeriodMatrix = new ComplexMatrix(new Complex[][] {{piNormalizedTau}});
				System.out.println("Tau Re " + tau.re);
				System.out.println("Tau Im " + tau.im);
				System.out.println("Tau Abs " + tau.abs());
				System.out.println("Tau Arg " + tau.arg());
				SiegelReduction siegel = new SiegelReduction(PeriodMatrix);
				System.out.println("After Siegel: " +  siegel.getReducedPeriodMatrix());
			} catch (Exception e) {
				e.printStackTrace();
			}
			fireJobProgress(1.0);
			break;
		// genus > 1 -------------------------
		default:
			System.out.println("unwrapping surface of genus " + genus + "...");
			if (usePetsc) {
				unwrapper = new HyperbolicUnwrapperPETSc();
			} else {
				unwrapper = new HyperbolicUnwrapper();
			}
			if (useSelectionCuts) unwrapper.setCutGraph(selectedEdges);
			if (!selectedVertices.isEmpty()) {
				CoVertex cutRoot = selectedVertices.iterator().next();
				unwrapper.setCutRoot(cutRoot);
			}
			unwrapper.setGradientTolerance(gradTolerance);
			unwrapper.setMaxIterations(maxIterations);
			unwrapper.unwrap(surface, genus, aSet);
			unwrapTime = System.currentTimeMillis();
			fireJobProgress(0.5);
			cutInfo = unwrapper.getCutInfo();
			lengthMap = unwrapper.getlengthMap();
			layoutRoot = unwrapper.getLayoutRoot();
			layoutTime = System.currentTimeMillis();
			fireJobProgress(1.0);
			break;
		}
		NumberFormat nf = new DecimalFormat("0.00");
		System.out.println("minimization took " + nf.format((unwrapTime - startTime) / 1000.0) + "sec.");
		System.out.println("layout took " + nf.format((layoutTime - unwrapTime) / 1000.0) + "sec.");
//		int brokenCount = 0;
//		int curveVertices = 0;
//		for (CoEdge e : surface.getEdges()) {
//			if (e.getAlpha() >= Math.PI) {
//				brokenCount++;
//			}
//		}
//		for (CoVertex v : surface.getVertices()) {
//			if (HalfEdgeUtils.isBoundaryVertex(v)) {
//				continue;
//			}
//			double sum = EuclideanLayout.calculateAngleSum(v);
//			if (Math.abs(sum - 2*Math.PI) > 1E-6) {
//				curveVertices++;
//			}
//		}
//		if (brokenCount > 0) {
//			System.err.println("WARNING: there are " + brokenCount + " broken triangles in the solution.");
//		}
//		if (curveVertices > 0) {
//			System.err.println("WARNING: there are " + curveVertices + " internal vertices with angle sum != 2PI.");
//		}
	}
	
	public CoHDS getSurface() {
		return surface;
	}

	public void setSurface(CoHDS surface) {
		this.surface = surface;
	}

	public boolean isUsePetsc() {
		return usePetsc;
	}

	public void setUsePetsc(boolean usePetsc) {
		this.usePetsc = usePetsc;
	}

	public int getNumCones() {
		return numCones;
	}

	public void setNumCones(int numCones) {
		this.numCones = numCones;
	}
	
	public void setQuantizationMode(QuantizationMode quantizationMode) {
		this.coneMode = quantizationMode;
	}
	
	public QuantizationMode getQuantizationMode() {
		return coneMode;
	}
	
	public void setToleranceExponent(double toleranceExp) {
		this.toleranceExp = toleranceExp;
	}
	
	public void setMaxIterations(int maxIterations) {
		this.maxIterations = maxIterations;
	}
	
	public void setBoundaryMode(BoundaryMode boundaryMode) {
		this.boundaryMode = boundaryMode;
	}
	
	public void setBoundaryQuantMode(QuantizationMode boundaryQuantMode) {
		this.boundaryQuantMode = boundaryQuantMode;
	}
	
	public void setSelectedVertices(Set<CoVertex> selectedVertices) {
		this.selectedVertices = selectedVertices;
	}
	
	public void setSelectedEdges(Set<CoEdge> selectedEdges) {
		this.selectedEdges = selectedEdges;
	}

	public void setSpherize(boolean spherize) {
		this.spherize = spherize;
	}
	public boolean isSpherize() {
		return spherize;
	}
	
	public void setUseSelectionCuts(boolean useSelectionCuts) {
		this.useSelectionCuts = useSelectionCuts;
	}
	public boolean isUseSelectionCuts() {
		return useSelectionCuts;
	}
	
}
