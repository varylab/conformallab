package de.varylab.discreteconformal.convergence;

import static de.varylab.discreteconformal.convergence.ConvergenceQuality.calculateQualityMeasure;
import static java.lang.Double.POSITIVE_INFINITY;
import static java.lang.Math.signum;

import java.util.HashSet;
import java.util.Set;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
import de.jreality.math.Pn;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.convergence.ConvergenceQuality.QualityMeasure;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.SphereUtility;
import de.varylab.discreteconformal.util.DiscreteEllipticUtility;

public class ConvergenceRandom extends ConvergenceSeries {

	private int 
		minExtraPoints = 0,
		maxExtraPoints = 0,
		incExtraPoints = 1,
		numOptSteps = 0;
	private boolean 	
		reuse = false,
		trianopt = false;
	private QualityMeasure
		qualityMeasure = QualityMeasure.MaxCrossRatio;
	private double
		qualityThreshold = POSITIVE_INFINITY,
		qualityExponent = 1.0;
	
	@Override
	protected OptionSet configureAndParseOptions(OptionParser p, String... args) {
		OptionSpec<String> qualityMeasureSpec = p.accepts("QM", "Mesh quality measure").withRequiredArg().defaultsTo("MaxCrossRatio");
		OptionSpec<Double> qualityExponentSpec = p.accepts("QE", "Quality measure exponent").withRequiredArg().ofType(Double.class).defaultsTo(1.0);
		OptionSpec<Double> qualityThresholdSpec = p.accepts("QT", "Quality measure threshold").withRequiredArg().ofType(Double.class).defaultsTo(POSITIVE_INFINITY);
		OptionSpec<Integer> minPointsSpec = p.accepts("min", "Minimum number of extra points").withRequiredArg().ofType(Integer.class).defaultsTo(0);
		OptionSpec<Integer> maxPointsSpec = p.accepts("max", "Maximum number of extra points").withRequiredArg().ofType(Integer.class).defaultsTo(50);
		OptionSpec<Integer> incPointsSpec = p.accepts("inc", "Increment").withRequiredArg().ofType(Integer.class).defaultsTo(1);
		OptionSpec<Integer> optStepsSpec = p.accepts("nopt", "Number of triangulation optimization steps").withRequiredArg().ofType(Integer.class).defaultsTo(1);
		p.accepts("reuse", "Reuse extra points");
		OptionSet opts = p.parse(args);
		
		minExtraPoints = minPointsSpec.value(opts);
		maxExtraPoints = maxPointsSpec.value(opts);
		incExtraPoints = incPointsSpec.value(opts);
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
		qualityExponent = qualityExponentSpec.value(opts);
		qualityThreshold = qualityThresholdSpec.value(opts);
		return opts;
	}
	
	@Override
	protected void perform() throws Exception {
		writeComment("index[1], absErr[2], argErr[3], reErr[4], imErr[5], quality[6], re[7], im[8]");
		for (int i = minExtraPoints; i <= maxExtraPoints; i += incExtraPoints) {
			if (reuse) rnd.setSeed(123);
			CoHDS hds = new CoHDS();
			// predefined vertices
			for (int vi = 0; vi < vertices.length; vi++) {
				CoVertex v = hds.addNewVertex();
				v.P = new double[] {vertices[vi][0], vertices[vi][1], vertices[vi][2], 1.0};
				Pn.setToLength(v.P, v.P, 1, Pn.EUCLIDEAN);
			}
			// additional points
			for (int j = 0; j < i; j++) {
				CoVertex v = hds.addNewVertex();
				v.P = new double[] {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian(), 1.0};
				Pn.setToLength(v.P, v.P, 1, Pn.EUCLIDEAN);
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
			double meshQuality = 0.0;
			try {
				Set<CoEdge> glueSet = new HashSet<CoEdge>();
				DiscreteEllipticUtility.generateEllipticImage(hds, 0, glueSet, branchIndices);
				meshQuality = calculateQualityMeasure(qualityMeasure, qualityExponent, hds);
				if (meshQuality > qualityThreshold) {
					System.out.println("Mesh discarded: Q=" + meshQuality + " > " + qualityThreshold);
					continue;
				}
//				HalfedgeIO.writeOBJ(hds, new AdapterSet(new PositionAdapter()), "data/tmp/tmp" + i + ".obj");
				tau = DiscreteEllipticUtility.calculateHalfPeriodRatio(hds, 1E-8);
			} catch (Exception e) {
				e.printStackTrace();
				continue;
			}
			
			// check if tau is on the wrong side of the fundamental domain
			double reErrCheck = Math.abs(tau.re - getExpectedTau().re);
			if (reErrCheck > 0.5) {
				tau.re -= signum(tau.re);
			}
			
			double absErr = tau.abs() - getExpectedTau().abs();
			double argErr = tau.arg() - getExpectedTau().arg();
			double reErr = tau.re - getExpectedTau().re;
			double imErr = tau.im - getExpectedTau().im;
			writeErrorLine(i + "\t" + absErr + "\t" + argErr + "\t" + reErr + "\t" + imErr + "\t" + meshQuality + "\t" + tau.re + "\t" + tau.im);
		}

	}

}
