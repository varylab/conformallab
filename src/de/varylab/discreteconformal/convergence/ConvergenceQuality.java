package de.varylab.discreteconformal.convergence;

import geom3d.Point;

import java.io.FileWriter;
import java.util.HashSet;
import java.util.Set;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
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
		trianopt = false;
	private QualityMeasure
		qualityMeasure = QualityMeasure.MaxCrossRatio;
	private double
		measureExponent = 1.0;
	
	
	public enum QualityMeasure {
		MaxCrossRatio,
		MeanCrossRatio,
		SumCrossRatio,
		MaxMultiRatio,
		MeanMultiRatio,
		SumMultiRatio,
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
		p.accepts("reuse", "Reuse extra points");
		
		OptionSet opts = p.parse(args);
		
		numSamples = numSamplesSpec.value(opts);
		numExtraPoints = numExtraPointsSpec.value(opts);
		reuse = opts.has("reuse");
		trianopt = opts.has("nopt");
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
		writeComment("Measure exponent: " + measureExponent);
		writeComment("index[1], absErr[2], argErr[3], reErr[4], imErr[5], qFun[6]");
		for (int i = 0; i < numSamples; i ++) {
			if (reuse) rnd.setSeed(123);
			CoHDS hds = new CoHDS();
			// predefined vertices
			for (int vi = 0; vi < vertices.length; vi++) {
				CoVertex v = hds.addNewVertex();
				v.getPosition().set(vertices[vi][0], vertices[vi][1], vertices[vi][2]);	
			}
			// extra points
			for (int j = 0; j < numExtraPoints; j++) {
				double[] pos = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
				Point pointPos = new Point(pos);
				pointPos.normalize();
				CoVertex v = hds.addNewVertex();
				v.setPosition(pointPos);
			}
			// optimize triangulation
			if (trianopt) {
				Set<CoVertex> fixedVertices = new HashSet<CoVertex>();
				fixedVertices.add(hds.getVertex(branchIndices[0]));
				fixedVertices.add(hds.getVertex(branchIndices[1]));
				fixedVertices.add(hds.getVertex(branchIndices[2]));
				fixedVertices.add(hds.getVertex(branchIndices[3]));
				SphereUtility.equalizeSphereVertices(hds, fixedVertices, numOptSteps, 1E-6);
			}
			Complex tau = null;
			try {
				Set<CoEdge> glueSet = new HashSet<CoEdge>();
				DiscreteEllipticUtility.generateEllipticImage(hds, 0, glueSet, branchIndices);
				tau = DiscreteEllipticUtility.calculateHalfPeriodRatio(hds, 1E-8);
			} catch (Exception e) {
				e.printStackTrace();
				continue;
			}
			double absErr = tau.abs() - getExpectedTau().abs();
			double argErr = tau.arg() - getExpectedTau().arg();
			double reErr = tau.re - getExpectedTau().re;
			double imErr = tau.im - getExpectedTau().im;
			double qFun = 0;
			
			switch (qualityMeasure) {
			case MaxCrossRatio:
				qFun = maxLengthCrossRatioFunction(hds, measureExponent);
				break;
			case MeanCrossRatio:
				qFun = meanLengthCrossRatioFunction(hds, measureExponent);
				break;
			case SumCrossRatio:
				qFun = sumLengthCrossRatioFunction(hds, measureExponent);
				break;
			case MaxMultiRatio:
				qFun = maxMultiRatioFunction(hds, measureExponent);
				break;
			case MeanMultiRatio:
				qFun = meanMultiRatioFunction(hds, measureExponent);
				break;
			case SumMultiRatio:
				qFun = sumMultiRatioFunction(hds, measureExponent);
				break;
			}
			
			writeErrorLine(i + "\t" + absErr + "\t" + argErr + "\t" + reErr + "\t" + imErr + "\t" + qFun);
		}
	}
	
	
	public double maxLengthCrossRatioFunction(CoHDS hds, double exp) {
		double r = 0.0;
		for (CoEdge e : hds.getPositiveEdges()) {
			double q = getLengthCrossRatio(e);
			double qfun = (q + 1/q)/2 - 1;
			r = Math.max(r, Math.pow(qfun, exp));
		}
		return r;
	}
	
	
	public double meanLengthCrossRatioFunction(CoHDS hds, double exp) {
		double r = 0.0;
		for (CoEdge e : hds.getPositiveEdges()) {
			double q = getLengthCrossRatio(e);
			double qfun = (q + 1/q)/2 - 1;
			r += Math.pow(qfun, exp);
		}
		return 2 * r / hds.numEdges();
	}
	
	public double sumLengthCrossRatioFunction(CoHDS hds, double exp) {
		double r = 0.0;
		for (CoEdge e : hds.getPositiveEdges()) {
			double q = getLengthCrossRatio(e);
			double qfun = (q + 1/q)/2 - 1;
			r += Math.pow(qfun, exp);
		}
		return r;
	}
	
	public double maxMultiRatioFunction(CoHDS hds, double exp) {
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
	
	public double meanMultiRatioFunction(CoHDS hds, double exp) {
		double r = 0.0;
		for (CoFace f : hds.getFaces()) {
			double q = getLengthMultiRatio(f);
			double qfun = (q + 1/q)/2 - 1;
			r += Math.pow(qfun, exp);
		}
		return r / hds.numFaces();
	}
	
	public double sumMultiRatioFunction(CoHDS hds, double exp) {
		double r = 0.0;
		for (CoFace f : hds.getFaces()) {
			double q = getLengthMultiRatio(f);
			double qfun = (q + 1/q)/2 - 1;
			r += Math.pow(qfun, exp);
		}
		return r;
	}
	
	
	private double getLengthCrossRatio(CoEdge e) {
		double a = e.getNextEdge().getLength();
		double b = e.getPreviousEdge().getLength();
		double c = e.getOppositeEdge().getNextEdge().getLength();
		double d = e.getOppositeEdge().getPreviousEdge().getLength();
		return (a * c) / (b * d);
	}
	
	private double getLengthMultiRatio(CoFace f) {
		double q = 1.0;
		for (CoEdge e : HalfEdgeUtils.boundaryEdges(f)) {
			q *= getLengthCrossRatio(e);
		}
		return q;
	}
	
	
}
