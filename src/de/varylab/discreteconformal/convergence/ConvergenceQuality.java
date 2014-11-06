package de.varylab.discreteconformal.convergence;

import static java.lang.Math.sqrt;

import java.io.FileWriter;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.SphereUtility;
import de.varylab.discreteconformal.util.DiscreteEllipticUtility;

public class ConvergenceQuality extends ConvergenceSeries {

	private int 
		numSamples = 1,
		numExtraPoints = 0,
		numOptSteps = 0;
	private boolean 
		reuse = false,
		trainoptrand = false,
		trianopt = false;
	private QualityMeasure
		qualityMeasure = QualityMeasure.MaxCrossRatio;
	private double
		measureExponent = 1.0;
	private Logger
		log = Logger.getLogger(ConvergenceQuality.class.getName());
	
	
	public static enum QualityMeasure {
		MaxCrossRatio,
		MeanCrossRatio,
		SumCrossRatio,
		MaxMultiRatio,
		MeanMultiRatio,
		SumMultiRatio,
		ElectrostaticEnergy,
		MaxCircumCircle
	}
	
	
	public ConvergenceQuality() {
	}


	public ConvergenceQuality(
		double[][] vertices, 
		int[] branchIndices,
		Complex tauExpected, 
		FileWriter errorWriter, 
		int numSamples,
		int numExtraPoints, 
		boolean reuse, 
		boolean trianopt,
		int numOptSteps
	) {
		super(vertices, branchIndices, tauExpected, errorWriter);
		this.numSamples = numSamples;
		this.numExtraPoints = numExtraPoints;
		this.numOptSteps = numOptSteps;
		this.reuse = reuse;
		this.trianopt = trianopt;
	}


	@Override
	protected OptionSet configureAndParseOptions(OptionParser p, String... args) {
		OptionSpec<String> qualityMeasureSpec = p.accepts("QM").withRequiredArg().defaultsTo("MaxCrossRatio");
		OptionSpec<Double> measureExpSpec = p.accepts("exp", "Measure function exponent").withRequiredArg().ofType(Double.class).defaultsTo(1.0);
		OptionSpec<Integer> numSamplesSpec = p.accepts("num", "Number of samples").withRequiredArg().ofType(Integer.class).defaultsTo(1);
		OptionSpec<Integer> numExtraPointsSpec = p.accepts("extra", "Number of extra points").withRequiredArg().ofType(Integer.class).defaultsTo(0);
		OptionSpec<Integer> optStepsSpec = p.accepts("nopt", "Number of triangulation optimization steps").withRequiredArg().ofType(Integer.class).defaultsTo(1);
		p.accepts("noptrand", "Randomize optimization steps");
		p.accepts("reuse", "Reuse extra points");
		
		OptionSet opts = p.parse(args);
		
		numSamples = numSamplesSpec.value(opts);
		numExtraPoints = numExtraPointsSpec.value(opts);
		reuse = opts.has("reuse");
		trianopt = opts.has("nopt");
		trainoptrand = opts.has("noptrand");
		numOptSteps = 0;
		if (trianopt) {
			numOptSteps = opts.valueOf(optStepsSpec);
			if (numOptSteps == 0) {
				trianopt = false;
			}
		}
		
		qualityMeasure = QualityMeasure.valueOf(qualityMeasureSpec.value(opts));
		if (qualityMeasure == null) {
			throw new RuntimeException("Unknown quality measure in ConvergenceQuality");
		}
		measureExponent = measureExpSpec.value(opts);
		
		return opts;
	}
	
	
	
	
	@Override
	protected void perform() throws Exception {
		writeComment("Quality measure: " + qualityMeasure);
		writeComment("Measure exponent: " + measureExponent);
		writeComment("index[1], absErr[2], argErr[3], reErr[4], imErr[5], re[6], im[7], qFun[8]");
		for (int i = 0; i < numSamples; i ++) {
			if (reuse) rnd.setSeed(123);
			CoHDS hds = new CoHDS();
			// predefined vertices
			for (int vi = 0; vi < vertices.length; vi++) {
				CoVertex v = hds.addNewVertex();
				v.P = new double[] {vertices[vi][0], vertices[vi][1], vertices[vi][2], 1.0};
				Pn.setToLength(v.P, v.P, 1.0, Pn.EUCLIDEAN);
			}
			// extra points
			for (int j = 0; j < numExtraPoints; j++) {
				CoVertex v = hds.addNewVertex();
				v.P = new double[] {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian(), 1.0};
				Pn.setToLength(v.P, v.P, 1.0, Pn.EUCLIDEAN);
			}
			// optimize triangulation
			if (trianopt) {
				Set<CoVertex> fixedVertices = new HashSet<CoVertex>();
				fixedVertices.add(hds.getVertex(branchIndices[0]));
				fixedVertices.add(hds.getVertex(branchIndices[1]));
				fixedVertices.add(hds.getVertex(branchIndices[2]));
				fixedVertices.add(hds.getVertex(branchIndices[3]));
				int steps = numOptSteps;
				if (trainoptrand) {
					steps = rnd.nextInt(numOptSteps + 1);
				}
				SphereUtility.equalizeSphereVertices(hds, fixedVertices, steps, 1E-6);
			}
			
			// calculate quality measure for the closed surface
			double meshQuality = 0;
			Complex tau = null;
			try {
				Set<CoEdge> glueSet = new HashSet<CoEdge>();
				DiscreteEllipticUtility.generateEllipticImage(hds, 0, glueSet, branchIndices);
				meshQuality = calculateQualityMeasure(qualityMeasure, measureExponent, hds);
				tau = DiscreteEllipticUtility.calculateHalfPeriodRatio(hds, 1E-9);
				if (qualityMeasure == QualityMeasure.MaxCircumCircle) {
					meshQuality = calculateQualityMeasure(qualityMeasure, measureExponent, hds);
				}
				if (Double.isInfinite(meshQuality) | Double.isNaN(meshQuality)) {
					continue;
				}
			} catch (Exception e) {
				e.printStackTrace();
				continue;
			}
			log.info("tau = " + tau);
			double absErr = tau.abs() - getExpectedTau().abs();
			double argErr = tau.arg() - getExpectedTau().arg();
			double reErr = tau.re - getExpectedTau().re;
			double imErr = tau.im - getExpectedTau().im;
			writeErrorLine(i + "\t" + absErr + "\t" + argErr + "\t" + reErr + "\t" + imErr + "\t" + tau.re + "\t" + tau.im + "\t" + meshQuality);
		}
	}
	
	
	
	public static double calculateQualityMeasure(QualityMeasure qualityMeasure, double measureExponent, CoHDS hds) {
		switch (qualityMeasure) {
		case MaxCrossRatio:
			return maxLengthCrossRatioFunction(hds, measureExponent);
		case MeanCrossRatio:
			return meanLengthCrossRatioFunction(hds, measureExponent);
		case SumCrossRatio:
			return sumLengthCrossRatioFunction(hds, measureExponent);
		case MaxMultiRatio:
			return maxMultiRatioFunction(hds, measureExponent);
		case MeanMultiRatio:
			return meanMultiRatioFunction(hds, measureExponent);
		case SumMultiRatio:
			return sumMultiRatioFunction(hds, measureExponent);
		case ElectrostaticEnergy:
			return electrostaticEnergy(hds);
		case MaxCircumCircle:
			return getMaxCircumRadius(hds) / getTextureArea(hds);
		}
		return 1.0;
	}
	
	
	public static double electrostaticEnergy(CoHDS hds) {
		double E = 0.0; 
		for (CoVertex v : hds.getVertices()) {
			double[] vPos = v.P;
			Pn.dehomogenize(vPos, vPos);
			for (CoVertex w : hds.getVertices()) {
				if (v == w) continue;
				double[] wPos = w.P;
				Pn.dehomogenize(wPos, wPos);
				double[] dir = Rn.subtract(null, vPos, wPos);
				double dsq = Rn.innerProduct(dir, dir);
				if (dsq == 0) continue;
				E += 1 / dsq;
			}
		}
		return E;
	}
	
	
	
	public static double maxLengthCrossRatioFunction(CoHDS hds, double exp) {
		double r = 0.0;
		for (CoEdge e : hds.getPositiveEdges()) {
			double q = getLengthCrossRatio(e);
			double qfun = (q + 1/q)/2 - 1;
			r = Math.max(r, Math.pow(qfun, exp));
		}
		return r;
	}
	
	
	public static double meanLengthCrossRatioFunction(CoHDS hds, double exp) {
		double r = 0.0;
		for (CoEdge e : hds.getPositiveEdges()) {
			double q = getLengthCrossRatio(e);
			double qfun = (q + 1/q)/2 - 1;
			r += Math.pow(qfun, exp);
		}
		return 2 * r / hds.numEdges();
	}
	
	public static double sumLengthCrossRatioFunction(CoHDS hds, double exp) {
		double r = 0.0;
		for (CoEdge e : hds.getPositiveEdges()) {
			double q = getLengthCrossRatio(e);
			double qfun = (q + 1/q)/2 - 1;
			r += Math.pow(qfun, exp);
		}
		return r;
	}
	
	public static double maxMultiRatioFunction(CoHDS hds, double exp) {
		double r = 0.0;
		for (CoFace f : hds.getFaces()) {
			double q = 1.0;
			for (CoEdge e : HalfEdgeUtils.boundaryEdges(f)) {
				q *= getLengthCrossRatio(e);
			}
			double qfun = (q + 1/q)/2 - 1;
			r = Math.max(r, Math.pow(qfun, exp));
		}
		return r;
	}
	
	public static double meanMultiRatioFunction(CoHDS hds, double exp) {
		double r = 0.0;
		for (CoFace f : hds.getFaces()) {
			double q = getLengthMultiRatio(f);
			double qfun = (q + 1/q)/2 - 1;
			r += Math.pow(qfun, exp);
		}
		return r / hds.numFaces();
	}
	
	public static double sumMultiRatioFunction(CoHDS hds, double exp) {
		double r = 0.0;
		for (CoFace f : hds.getFaces()) {
			double q = getLengthMultiRatio(f);
			double qfun = (q + 1/q)/2 - 1;
			r += Math.pow(qfun, exp);
		}
		return r;
	}
	
	public static double getMaxCircumRadius(CoHDS hds) {
		double r = 0.0;
		for (CoFace f : hds.getFaces()) {
			double rad = getTextureCircumCircleRadius(f);
			if (rad > r) {
				r = rad;
			}
		}
		return r;
	}
	
	public static double getTextureCircumCircleRadius(CoFace f) {
		CoEdge e = f.getBoundaryEdge();
		double a = e.getTexLength();
		double b = e.getNextEdge().getTexLength();
		double c = e.getPreviousEdge().getTexLength();
		double A = getTextureTriangleArea(f);
		return a*b*c / A / 4;
	}
	
	public static double getTextureArea(CoHDS hds) {
		double A = 0.0;
		for (CoFace f : hds.getFaces()) {
			A += getTextureTriangleArea(f);
		}
		return A;
	}
	
	public static double getTextureTriangleArea(CoFace f) {
		CoEdge e = f.getBoundaryEdge();
		double a = e.getTexLength();
		double b = e.getNextEdge().getTexLength();
		double c = e.getPreviousEdge().getTexLength();
		return sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c)) / 4;	
	}
	
	
	private static double getLengthCrossRatio(CoEdge e) {
		double a = e.getNextEdge().getLength();
		double b = e.getPreviousEdge().getLength();
		double c = e.getOppositeEdge().getNextEdge().getLength();
		double d = e.getOppositeEdge().getPreviousEdge().getLength();
		return (a * c) / (b * d);
	}
	
	private static double getLengthMultiRatio(CoFace f) {
		double q = 1.0;
		for (CoEdge e : HalfEdgeUtils.boundaryEdges(f)) {
			q *= getLengthCrossRatio(e);
		}
		return q;
	}
	
	
}
