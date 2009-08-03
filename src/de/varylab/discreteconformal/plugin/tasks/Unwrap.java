package de.varylab.discreteconformal.plugin.tasks;

import static de.varylab.discreteconformal.util.CuttingUtility.cutManifoldToDisk;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Collection;

import javax.swing.SwingWorker;

import no.uib.cipr.matrix.Vector;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.adapter.EuclideanLengthWeightAdapter;
import de.varylab.discreteconformal.adapter.HyperbolicLengthWeightAdapter;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.CHyperbolicLayout;
import de.varylab.discreteconformal.unwrapper.CHyperbolicUnwrapper;
import de.varylab.discreteconformal.unwrapper.CHyperbolicUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.EuclideanLayout;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapper;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.Unwrapper;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class Unwrap extends SwingWorker<CoHDS, Object> {

	private CoHDS
		surface = null;
	private boolean
		usePetsc = false,
		quantizeCones = false;
	private QuantizationMode
		quantizationMode = QuantizationMode.Quads;
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
	
	
	public static enum QuantizationMode {
		Quads,
		Hexagons
	}
	
	
	public Unwrap(CoHDS surface) {
		this.surface = surface;
	}
	
	
	@Override
	protected CoHDS doInBackground() throws Exception {
		long startTime = System.currentTimeMillis();
		long unwrapTime = -1;
		long layoutTime = -1;
		setProgress(0);
		genus = HalfEdgeUtils.getGenus(surface);
		Vector u = null;
		double gradTolerance = Math.pow(10, toleranceExp);
		switch (genus) {
		// disk or sphere---------------------
		case 0: 
			Collection<CoEdge> bList = HalfEdgeUtils.boundaryEdges(surface);
			boolean isSphere = bList.size() == 0;
			if (isSphere) {
				System.out.println("unwrapping a sphere...");
			} else {
				System.out.println("unwrapping a disk...");
			}
			Unwrapper unwrapper = null;
			if (usePetsc) {
				EuclideanUnwrapperPETSc uw = new EuclideanUnwrapperPETSc();
				uw.setNumCones(numCones);
				uw.setQuantizeCones(quantizeCones);
				uw.setQuantizationMode(quantizationMode);
				unwrapper = uw;
			} else {
				EuclideanUnwrapper uw = new EuclideanUnwrapper();
				uw.setNumCones(numCones);
				uw.setQuantizeCones(quantizeCones);
				uw.setQuantizationMode(quantizationMode);
				unwrapper = uw;
			}
			unwrapper.setGradientTolerance(gradTolerance);
			unwrapper.setMaxIterations(maxIterations);
			u = unwrapper.unwrap(surface);
			System.err.println("---minimzation redady---");
			unwrapTime = System.currentTimeMillis();
			setProgress(50);
			layoutRoot = EuclideanLayout.doLayout(surface, u);
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
			u = unwrapper.unwrap(surface);
			unwrapTime = System.currentTimeMillis();
			setProgress(50);
			EuclideanLengthWeightAdapter eucWa = new EuclideanLengthWeightAdapter(u); 
			CoVertex cutRoot = surface.getVertex(getMinUIndex(u));
			cutInfo = cutManifoldToDisk(surface, cutRoot, eucWa);
			layoutRoot = EuclideanLayout.doLayout(surface, u);
			layoutTime = System.currentTimeMillis();
			setProgress(100);
			break;
		// genus > 1 -------------------------
		default:
			System.out.println("unwrapping surface of genus " + genus + "...");
			if (usePetsc) {
				unwrapper = new CHyperbolicUnwrapperPETSc();
			} else {
				unwrapper = new CHyperbolicUnwrapper();
			}
			unwrapper.setGradientTolerance(gradTolerance);
			unwrapper.setMaxIterations(maxIterations);
			u = unwrapper.unwrap(surface);
			System.err.println("---minimization ready---");
			unwrapTime = System.currentTimeMillis();
			setProgress(50);
			HyperbolicLengthWeightAdapter hypWa = new HyperbolicLengthWeightAdapter(u);
			cutRoot = surface.getVertex(getMinUIndex(u));
			cutInfo = cutManifoldToDisk(surface, cutRoot, hypWa);
			CoVertex layoutRoot = surface.getVertex(getMaxUIndex(u));
			layoutRoot = CHyperbolicLayout.doLayout(surface, layoutRoot, u);
			layoutTime = System.currentTimeMillis();
			setProgress(100);
			break;
		}
		NumberFormat nf = new DecimalFormat("0.00");
		System.out.println("minimization took " + nf.format((unwrapTime - startTime) / 1000.0) + "sec.");
		System.out.println("layout took " + nf.format((layoutTime - unwrapTime) / 1000.0) + "sec.");
		return surface;
	}
	
	
	private int getMinUIndex(Vector u) {
		int index = 0;
		double iVal = u.get(0);
		for (int i = 1; i < u.size(); i++) {
			double val = u.get(i);
			if (iVal < val) {
				index = i;
				iVal = i;
			}
		}
		return index;
	}
	
	private int getMaxUIndex(Vector u) {
		int index = 0;
		double iVal = u.get(0);
		for (int i = 1; i < u.size(); i++) {
			double val = u.get(i);
			if (iVal > val) {
				index = i;
				iVal = i;
			}
		}
		return index;
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


	public boolean isQuantizeCones() {
		return quantizeCones;
	}


	public void setQuantizeCones(boolean quantizeCones) {
		this.quantizeCones = quantizeCones;
	}


	public int getNumCones() {
		return numCones;
	}


	public void setNumCones(int numCones) {
		this.numCones = numCones;
	}
	
	public void setQuantizationMode(QuantizationMode quantizationMode) {
		this.quantizationMode = quantizationMode;
	}
	
	public QuantizationMode getQuantizationMode() {
		return quantizationMode;
	}
	
	public void setToleranceExponent(double toleranceExp) {
		this.toleranceExp = toleranceExp;
	}
	
	public void setMaxIterations(int maxIterations) {
		this.maxIterations = maxIterations;
	}

}
