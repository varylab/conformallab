package de.varylab.discreteconformal.convergence;

import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
import de.jreality.math.Pn;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.SphereUtility;
import de.varylab.discreteconformal.util.DiscreteEllipticUtility;

public class ConvergenceNumberOfPoints extends ConvergenceSeries {

	private Logger	
		log = Logger.getLogger(ConvergenceNumberOfPoints.class.getName());
	private int 
		minExtraPoints = 0,
		maxExtraPoints = 0,
		iterations = 1,
		numOptSteps = 0,
		numJobs = 4;
	
	@Override
	protected OptionSet configureAndParseOptions(OptionParser p, String... args) {
		OptionSpec<Integer> minPointsSpec = p.accepts("min", "Minimum number of extra points").withRequiredArg().ofType(Integer.class).defaultsTo(0);
		OptionSpec<Integer> maxPointsSpec = p.accepts("max", "Maximum number of extra points").withRequiredArg().ofType(Integer.class).defaultsTo(50);
		OptionSpec<Integer> iterationsPointsSpec = p.accepts("iter", "Iterations").withRequiredArg().ofType(Integer.class).defaultsTo(1000);
		OptionSpec<Integer> optStepsSpec = p.accepts("nopt", "Number of triangulation optimization steps").withRequiredArg().ofType(Integer.class).defaultsTo(0);
		OptionSpec<Integer> jobsSpec = p.accepts("jobs", "Number of Jobs").withRequiredArg().ofType(Integer.class).defaultsTo(4);
		OptionSet opts = p.parse(args);
		
		minExtraPoints = minPointsSpec.value(opts);
		maxExtraPoints = maxPointsSpec.value(opts);
		iterations = iterationsPointsSpec.value(opts);
		numOptSteps = opts.valueOf(optStepsSpec);
		numJobs = opts.valueOf(jobsSpec);
		return opts;
	}
	
	@Override
	protected void perform() throws Exception {
		String description = "vertexNum[1]\tdistErr[2]\tabsErr[3]\targErr[4]\treErr[5]\timErr[6]\tre[7]\tim[8]\t"; 
		description += "MaxCrossRatio\t";
		description += "MeanCrossRatio\t";
		description += "SumCrossRatio\t";
		description += "MaxMultiRatio\t";
		description += "MeanMultiRatio\t";
		description += "SumMultiRatio\t";
		description += "MaxCircleRadius\t";
		description += "MeanCircleRadius\t";
		description += "SumCircleRadius";
		writeComment(description);
		ExecutorService executor = Executors.newFixedThreadPool(numJobs);
		for (int i = 0; i < iterations; i++) {
			executor.submit(new Runnable() {
				@Override
				public void run() {
					CoHDS hds = new CoHDS();
					// predefined vertices
					for (int vi = 0; vi < vertices.length; vi++) {
						CoVertex v = hds.addNewVertex();
						v.P = new double[] {vertices[vi][0], vertices[vi][1], vertices[vi][2], 1.0};
						Pn.setToLength(v.P, v.P, 1, Pn.EUCLIDEAN);
					}
					// additional points
					int numVertices = rnd.nextInt(maxExtraPoints - minExtraPoints) + minExtraPoints;
					for (int j = 0; j < numVertices; j++) {
						CoVertex v = hds.addNewVertex();
						v.P = new double[] {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian(), 1.0};
						Pn.setToLength(v.P, v.P, 1, Pn.EUCLIDEAN);
					}
					// optimize triangulation
					if (numOptSteps > 0) {
						Set<CoVertex> fixedVertices = new HashSet<CoVertex>();
						fixedVertices.add(hds.getVertex(branchIndices[0]));
						fixedVertices.add(hds.getVertex(branchIndices[1]));
						fixedVertices.add(hds.getVertex(branchIndices[2]));
						fixedVertices.add(hds.getVertex(branchIndices[3]));
						SphereUtility.equalizeSphereVertices(hds, fixedVertices, numOptSteps, 1E-6);
					}
					Complex tau = null;
					double[] crossRatioQuality = null;
					double[] multiRatioQuality = null;
					double[] circleRadiusQuality = null;
					try {
						Set<CoEdge> glueSet = new HashSet<CoEdge>();
						DiscreteEllipticUtility.generateEllipticImage(hds, 0, glueSet, branchIndices);
						crossRatioQuality = ConvergenceUtility.getMaxMeanSumCrossRatio(hds, 1);
						multiRatioQuality = ConvergenceUtility.getMaxMeanSumMultiRatio(hds, 1);
						tau = DiscreteEllipticUtility.calculateHalfPeriodRatio(hds, 1E-9);
						circleRadiusQuality = ConvergenceUtility.getMaxMeanSumScaleInvariantCircumRadius(hds);
					} catch (Exception e) {
						log.warning(e.getMessage());
						return;
					}
					
					double absErr = tau.abs() - getExpectedTau().abs();
					double argErr = tau.arg() - getExpectedTau().arg();
					double reErr = tau.re - getExpectedTau().re;
					double imErr = tau.im - getExpectedTau().im;
					double distErr = getExpectedTau().minus(tau).abs();
					String logString = numVertices + "\t" + distErr + "\t" + absErr + "\t" + argErr + "\t" + reErr + "\t" + imErr + "\t" + tau.re + "\t" + tau.im + "\t";
					logString += crossRatioQuality[0] + "\t";
					logString += crossRatioQuality[1] + "\t";
					logString += crossRatioQuality[2] + "\t";
					logString += multiRatioQuality[0] + "\t";
					logString += multiRatioQuality[1] + "\t";
					logString += multiRatioQuality[2] + "\t";
					logString += circleRadiusQuality[0] + "\t";
					logString += circleRadiusQuality[1] + "\t";
					logString += circleRadiusQuality[2] + "\t";
					writeData(logString);
				}
			});
		}
		executor.awaitTermination(10, TimeUnit.DAYS);
	}

}
