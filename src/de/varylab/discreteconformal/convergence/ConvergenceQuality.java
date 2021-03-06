package de.varylab.discreteconformal.convergence;

import static de.varylab.discreteconformal.convergence.ConvergenceUtility.getMaxMeanSumScaleInvariantCircumRadius;
import static de.varylab.discreteconformal.convergence.ConvergenceUtility.getMaxMeanSumCrossRatio;
import static de.varylab.discreteconformal.convergence.ConvergenceUtility.getMaxMeanSumMultiRatio;
import static de.varylab.discreteconformal.util.DiscreteEllipticUtility.calculateHalfPeriodRatio;
import static de.varylab.discreteconformal.util.DiscreteEllipticUtility.generateEllipticImage;

import java.io.FileWriter;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
import de.jreality.math.Pn;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.SphereUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class ConvergenceQuality extends ConvergenceSeries {

	private int 
		numSamples = 1,
		minExtraPoints = 0,
		numExtraPoints = 0,
		numOptSteps = 0;
	private boolean 
		numExtraRand = false,
		reuse = false,
		trianoptrand = false,
		trianopt = false;
	private double
		measureExponent = 1.0;
	private Logger
		log = Logger.getLogger(ConvergenceQuality.class.getName());
	
	
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
		OptionSpec<Double> measureExpSpec = p.accepts("exp", "Measure function exponent").withRequiredArg().ofType(Double.class).defaultsTo(1.0);
		OptionSpec<Integer> numSamplesSpec = p.accepts("num", "Number of samples").withRequiredArg().ofType(Integer.class).defaultsTo(1);
		OptionSpec<Integer> numExtraPointsSpec = p.accepts("extra", "Number of extra points").withRequiredArg().ofType(Integer.class).defaultsTo(0);
		OptionSpec<Integer> minExtraPointsSpec = p.accepts("minextra", "Minimum number of extra points when randomizing").withRequiredArg().ofType(Integer.class).defaultsTo(0);
		OptionSpec<Integer> optStepsSpec = p.accepts("nopt", "Number of triangulation optimization steps").withRequiredArg().ofType(Integer.class).defaultsTo(1);
		p.accepts("noptrand", "Randomize optimization steps");
		p.accepts("extrarand", "Randomize number of extra points");
		p.accepts("reuse", "Reuse extra points");
		
		OptionSet opts = p.parse(args);
		
		numSamples = numSamplesSpec.value(opts);
		numExtraPoints = numExtraPointsSpec.value(opts);
		reuse = opts.has("reuse");
		trianopt = opts.has("nopt");
		trianoptrand = opts.has("noptrand");
		numExtraRand = opts.has("extrarand");
		minExtraPoints = minExtraPointsSpec.value(opts);
		numOptSteps = 0;
		if (trianopt) {
			numOptSteps = opts.valueOf(optStepsSpec);
			if (numOptSteps == 0) {
				trianopt = false;
			}
		}
		measureExponent = measureExpSpec.value(opts);
		return opts;
	}
	
	
	
	
	@Override
	protected void perform() throws Exception {
		writeComment("Measure exponent: " + measureExponent);
		String description = "index[1]\tnumVertices[2]\tdistErr[3]\tabsErr[4]\targErr[5]\treErr[6]\timErr[7]\tre[8]\tim[9]\t";
		description += "MaxCrossRatio[10]\t";
		description += "MeanCrossRatio[11]\t";
		description += "SumCrossRatio[12]\t";
		description += "MaxMultiRatio[13]\t";
		description += "MeanMultiRatio[14]\t";
		description += "SumMultiRatio[15]\t";
		description += "MaxCircleRadius[16]\t";
		description += "MeanCircleRadius[17]\t";
		description += "SumCircleRadius[18]";
		writeComment(description);
		for (int i = 0; i < numSamples; i ++) {
			if (reuse) rnd.setSeed(123);
			CoHDS hds = new CoHDS();
			// predefined vertices
			for (int vi = 0; vi < vertices.length; vi++) {
				CoVertex v = hds.addNewVertex();
				v.P = new double[] {vertices[vi][0], vertices[vi][1], vertices[vi][2], 1.0};
				Pn.setToLength(v.P, v.P, 1.0, Pn.EUCLIDEAN);
			}
			// extra points, possibly a random number of points
			int numExtra = numExtraPoints;
			if (numExtraRand) {
				numExtra = rnd.nextInt(numExtraPoints - minExtraPoints) + minExtraPoints;
			}
			for (int j = 0; j < numExtra; j++) {
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
				if (trianoptrand) {
					steps = rnd.nextInt(numOptSteps + 1);
				}
				SphereUtility.equalizeSphereVertices(hds, fixedVertices, steps, 1E-6);
			}
			
			// calculate quality measure for the closed surface
			Complex tau = null;
			double[] crossRatioQuality = null;
			double[] multiRatioQuality = null;
			double[] circleRadiusQuality = null;
			CoVertex cutRoot = hds.getVertex(branchIndices[0]);
			CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<>();
			int numVertices = 0;
			try {
				Set<CoEdge> glueSet = new HashSet<CoEdge>();
				Map<CoVertex, CoVertex> involution = generateEllipticImage(hds, 0, glueSet, branchIndices);
				if (!cutRoot.isValid()) cutRoot = involution.get(cutRoot);
				numVertices = hds.numVertices();
				crossRatioQuality = getMaxMeanSumCrossRatio(hds, measureExponent);
				multiRatioQuality = getMaxMeanSumMultiRatio(hds, measureExponent);
				tau = calculateHalfPeriodRatio(hds, cutRoot, 1E-9, cutInfo);
				circleRadiusQuality = getMaxMeanSumScaleInvariantCircumRadius(hds);
			} catch (Exception e) {
				e.printStackTrace();
				continue;
			}
			log.info("tau = " + tau);
			double absErr = tau.abs() - getExpectedTau().abs();
			double argErr = tau.arg() - getExpectedTau().arg();
			double reErr = tau.re - getExpectedTau().re;
			double imErr = tau.im - getExpectedTau().im;
			double distErr = getExpectedTau().minus(tau).abs();
			String logString = i + "\t" + numVertices + "\t" + distErr + "\t" + absErr + "\t" + argErr + "\t" + reErr + "\t" + imErr + "\t" + tau.re + "\t" + tau.im + "\t";
			logString += crossRatioQuality[0] + "\t";
			logString += crossRatioQuality[1] + "\t";
			logString += crossRatioQuality[2] + "\t";
			logString += multiRatioQuality[0] + "\t";
			logString += multiRatioQuality[1] + "\t";
			logString += multiRatioQuality[2] + "\t";
			logString += circleRadiusQuality[0] + "\t";
			logString += circleRadiusQuality[1] + "\t";
			logString += circleRadiusQuality[2];
			writeData(logString);
		}
	}
	
}
