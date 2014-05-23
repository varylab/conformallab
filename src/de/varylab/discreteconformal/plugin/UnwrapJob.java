package de.varylab.discreteconformal.plugin;

import static de.varylab.discreteconformal.plugin.DiscreteConformalPlugin.CHANNEL_BROKEN_TRIANGLES;
import static java.lang.Math.PI;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

import de.jreality.plugin.job.AbstractJob;
import de.jtem.blas.ComplexMatrix;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Selection;
import de.jtem.halfedgetools.selection.TypedSelection;
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
import de.varylab.discreteconformal.util.UnwrapUtility;

public class UnwrapJob extends AbstractJob {

	private static Logger 
		log = Logger.getLogger(UnwrapJob.class.getName());
	
	private CoHDS
		surface = null;
	private AdapterSet
		aSet = new AdapterSet();
	private TypedSelection<CoVertex>
		selectedVertices = new TypedSelection<CoVertex>();
	private TypedSelection<CoEdge>
		selectedEdges = new TypedSelection<CoEdge>();
	private TargetGeometry
		targetGeometry = TargetGeometry.Automatic;
	private boolean
		usePetsc = false,
		useSelectionCuts = false;
	private QuantizationMode
		coneMode = QuantizationMode.AllAngles,
		boundaryQuantizationMode = QuantizationMode.AllAngles;
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
	
	
	
	public UnwrapJob(CoHDS surface, AdapterSet aSet) {
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
		double gradTolerance = Math.pow(10, toleranceExp);
		
		if (targetGeometry == TargetGeometry.Automatic) {
			targetGeometry = UnwrapUtility.calculateTargetGeometry(surface);
		}
		
		switch (targetGeometry) {
		case Spherical: 
			log.info("Spherical unwrap...");
			Unwrapper unwrapper = null;
			if (usePetsc) {
//				unwrapper = new SphericalUnwrapperPETSc();
				unwrapper = new StereographicUnwrapper();
			} else {
//				unwrapper = new SphericalUnwrapper();
				unwrapper = new StereographicUnwrapper();
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
					log.warning("face " + f + " has an area greater than 2PI");
				}
			}
			cutInfo = unwrapper.getCutInfo();
			lengthMap = unwrapper.getlengthMap();
			layoutRoot = unwrapper.getLayoutRoot();
			layoutTime = System.currentTimeMillis();
			fireJobProgress(1.0);
			break;
		case Euclidean:
			log.info("Euclidean unwrap...");
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
					uw.setBoundaryQuantMode(boundaryQuantizationMode);
					unwrapper = uw;
				} else {
					EuclideanUnwrapper uw = new EuclideanUnwrapper();
					uw.setNumCones(numCones);
					uw.setConeMode(coneMode);
					uw.setBoundaryMode(boundaryMode);
					uw.setBoundaryQuantMode(boundaryQuantizationMode);
					unwrapper = uw;
				}
			}
			if (useSelectionCuts) {
				Set<CoEdge> edgeSet = new TreeSet<>(selectedEdges);
				unwrapper.setCutGraph(edgeSet);
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
			
			cutInfo = unwrapper.getCutInfo();
			lengthMap = unwrapper.getlengthMap();
			layoutRoot = unwrapper.getLayoutRoot();
			layoutTime = System.currentTimeMillis();
			if (genus == 1) {
				try {
					Complex tau = DiscreteEllipticUtility.calculateHalfPeriodRatio(cutInfo);
					Complex piNormalizedTau = new Complex(0, 2 * PI);
					piNormalizedTau = piNormalizedTau.times(tau);
					final ComplexMatrix PeriodMatrix = new ComplexMatrix(new Complex[][] {{piNormalizedTau}});
					log.info("Tau Re " + tau.re);
					log.info("Tau Im " + tau.im);
					log.info("Tau Abs " + tau.abs());
					log.info("Tau Arg " + tau.arg());
					SiegelReduction siegel = new SiegelReduction(PeriodMatrix);
					log.info("After Siegel: " +  siegel.getReducedPeriodMatrix());
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			fireJobProgress(1.0);
			break;
		default:
			log.info("Hyperbolic unwrap...");
			if (usePetsc) {
				unwrapper = new HyperbolicUnwrapperPETSc();
			} else {
				unwrapper = new HyperbolicUnwrapper();
			}
			if (useSelectionCuts) {
				Set<CoEdge> edgeSet = new TreeSet<>(selectedEdges);
				unwrapper.setCutGraph(edgeSet);
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
			cutInfo = unwrapper.getCutInfo();
			lengthMap = unwrapper.getlengthMap();
			layoutRoot = unwrapper.getLayoutRoot();
			layoutTime = System.currentTimeMillis();
			fireJobProgress(1.0);
			break;
		}
		NumberFormat nf = new DecimalFormat("0.00");
		log.info("minimization took " + nf.format((unwrapTime - startTime) / 1000.0) + "sec.");
		log.info("layout took " + nf.format((layoutTime - unwrapTime) / 1000.0) + "sec.");
		int brokenCount = 0;
		int curveVertices = 0;
		for (CoEdge e : surface.getEdges()) {
			if (e.getLeftFace() == null) continue;
			if (e.getAlpha() >= Math.PI) {
				brokenCount++;
				aSet.set(Selection.class, e.getLeftFace(), CHANNEL_BROKEN_TRIANGLES);
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
			log.warning("there are " + brokenCount + " broken triangles in the solution.");
		}
		if (curveVertices > 0) {
			log.warning("there are " + curveVertices + " internal vertices with angle sum != 2PI.");
		}
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
	
	public BoundaryMode getBoundaryMode() {
		return boundaryMode;
	}
	public void setBoundaryMode(BoundaryMode boundaryMode) {
		this.boundaryMode = boundaryMode;
	}
	
	public void setBoundaryQuantizationMode(QuantizationMode boundaryQuantMode) {
		this.boundaryQuantizationMode = boundaryQuantMode;
	}
	
	public void setSelectedVertices(TypedSelection<CoVertex> selectedVertices) {
		this.selectedVertices = selectedVertices;
	}
	
	public void setSelectedEdges(TypedSelection<CoEdge> selectedEdges) {
		this.selectedEdges = selectedEdges;
	}

	public void setUseSelectionCuts(boolean useSelectionCuts) {
		this.useSelectionCuts = useSelectionCuts;
	}
	public boolean isUseSelectionCuts() {
		return useSelectionCuts;
	}
	
	public void setTargetGeometry(TargetGeometry targetGeometry) {
		this.targetGeometry = targetGeometry;
	}
	public TargetGeometry getTargetGeometry() {
		return targetGeometry;
	}
	
}
