package de.varylab.discreteconformal.plugin.tasks;

import static java.lang.Math.PI;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.swing.SwingWorker;

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
import de.varylab.discreteconformal.unwrapper.EuclideanLayout;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapper;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.HyperbolicUnwrapper;
import de.varylab.discreteconformal.unwrapper.HyperbolicUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.QuantizationMode;
import de.varylab.discreteconformal.unwrapper.StereographicUnwrapper;
import de.varylab.discreteconformal.unwrapper.Unwrapper;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.DiscreteEllipticUtility;

public class Unwrap extends SwingWorker<CoHDS, Void> {

	private CoHDS
		surface = null;
	private AdapterSet
		aSet = new AdapterSet();
	private Set<CoVertex>
		selectedVertices = new HashSet<CoVertex>();
	private boolean
		usePetsc = false;
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
	protected CoHDS doInBackground() throws Exception {
		long startTime = System.currentTimeMillis();
		long unwrapTime = -1;
		long layoutTime = -1;
		setProgress(0);
		if (!HalfEdgeUtils.isValidSurface(surface, true)) {
			throw new RuntimeException("Surface is not valid");
		}
		genus = HalfEdgeUtils.getGenus(surface);
		double gradTolerance = Math.pow(10, toleranceExp);
		switch (genus) {
		// disk or sphere---------------------
		case 0: 
			Collection<CoEdge> bList = HalfEdgeUtils.boundaryEdges(surface);
			boolean isSphere = bList.size() == 0;
			Unwrapper unwrapper = null;
			if (isSphere) {
				System.out.println("unwrapping a sphere...");
				unwrapper = new StereographicUnwrapper();
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
			setProgress(50);
			
			cutInfo = unwrapper.getCutInfo();
			lengthMap = unwrapper.getlengthMap();
			layoutRoot = unwrapper.getLayoutRoot();
			layoutTime = System.currentTimeMillis();
			setProgress(100);
			break;
		// torus ----------------------------
		case 1:
			System.out.println("unwrapping a torus...");
			if (usePetsc) {
				unwrapper = new EuclideanUnwrapperPETSc();
			} else {
				unwrapper = new EuclideanUnwrapper();
			}
			unwrapper.setGradientTolerance(gradTolerance);
			unwrapper.setMaxIterations(maxIterations);
			unwrapper.unwrap(surface, genus, aSet);
			unwrapTime = System.currentTimeMillis();
			setProgress(50);
			
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
			setProgress(100);
			break;
		// genus > 1 -------------------------
		default:
			System.out.println("unwrapping surface of genus " + genus + "...");
			if (usePetsc) {
				HyperbolicUnwrapperPETSc uw = new HyperbolicUnwrapperPETSc();
				if (!selectedVertices.isEmpty()) {
					CoVertex cutRoot = selectedVertices.iterator().next();
					uw.setCutRoot(cutRoot);
				}
				unwrapper = uw;
			} else {
				HyperbolicUnwrapper uw = new HyperbolicUnwrapper();
				if (!selectedVertices.isEmpty()) {
					CoVertex cutRoot = selectedVertices.iterator().next();
					uw.setCutRoot(cutRoot);
				}
				unwrapper = uw;
			}
			unwrapper.setGradientTolerance(gradTolerance);
			unwrapper.setMaxIterations(maxIterations);
			unwrapper.unwrap(surface, genus, aSet);
			unwrapTime = System.currentTimeMillis();
			setProgress(50);
			cutInfo = unwrapper.getCutInfo();
			lengthMap = unwrapper.getlengthMap();
			layoutRoot = unwrapper.getLayoutRoot();
			layoutTime = System.currentTimeMillis();
			setProgress(100);
			break;
		}
		NumberFormat nf = new DecimalFormat("0.00");
		System.out.println("minimization took " + nf.format((unwrapTime - startTime) / 1000.0) + "sec.");
		System.out.println("layout took " + nf.format((layoutTime - unwrapTime) / 1000.0) + "sec.");
		int brokenCount = 0;
		int curveVertices = 0;
		for (CoEdge e : surface.getEdges()) {
			if (e.getAlpha() >= Math.PI) {
				brokenCount++;
			}
		}
		for (CoVertex v : surface.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) {
				continue;
			}
			double sum = EuclideanLayout.calculateAngleSum(v);
			if (Math.abs(sum - 2*Math.PI) > 1E-6) {
				curveVertices++;
			}
		}
		if (brokenCount > 0) {
			System.err.println("WARNING: there are " + brokenCount + " broken triangles in the solution.");
		}
		if (curveVertices > 0) {
			System.err.println("WARNING: there are " + curveVertices + " internal vertices with angle sum != 2PI.");
		}
		return surface;
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

}
