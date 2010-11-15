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
import de.varylab.discreteconformal.util.DiscreteEllipticUtility;

public class ConvergenceNoise extends ConvergenceSeries {

	private QualityMeasure
		qualityMeasure = QualityMeasure.MaxCrossRatio;
	private double
		qualityThreshold = POSITIVE_INFINITY,
		qualityExponent = 1.0,
		noiseCoeff = 1E-8;
	private int
		numIterations = 1;
	
	@Override
	protected OptionSet configureAndParseOptions(OptionParser p, String... args) {
		OptionSpec<String> qualityMeasureSpec = p.accepts("QM", "Mesh quality measure").withRequiredArg().defaultsTo("MaxCrossRatio");
		OptionSpec<Double> qualityExponentSpec = p.accepts("QE", "Quality measure exponent").withRequiredArg().ofType(Double.class).defaultsTo(1.0);
		OptionSpec<Double> qualityThresholdSpec = p.accepts("QT", "Quality measure threshold").withRequiredArg().ofType(Double.class).defaultsTo(POSITIVE_INFINITY);
		OptionSpec<Integer> numIterationsSpec = p.accepts("IT", "Number of iterations").withRequiredArg().ofType(Integer.class).defaultsTo(1);
		OptionSpec<Double> noiseCoefficientSpec = p.accepts("NC", "Noise coefficient").withRequiredArg().ofType(Double.class).defaultsTo(1E-8);
		OptionSet opts = p.parse(args);
		
		numIterations = numIterationsSpec.value(opts);
		noiseCoeff = noiseCoefficientSpec.value(opts);
		qualityMeasure = QualityMeasure.valueOf(qualityMeasureSpec.value(opts));
		qualityExponent = qualityExponentSpec.value(opts);
		qualityThreshold = qualityThresholdSpec.value(opts);
		return opts;
	}
	
	@Override
	protected void perform() throws Exception {
		writeComment("index[1], absErr[2], argErr[3], reErr[4], imErr[5], quality[6]");
		
		double[] vPos = {0,0,0,1};
		for (int i = 0; i < numIterations; i ++) {
			CoHDS hds = new CoHDS();
			// predefined vertices
			for (int vi = 0; vi < vertices.length; vi++) {
				CoVertex v = hds.addNewVertex();
				vPos[0] = vertices[vi][0] + noiseCoeff * rnd.nextDouble();
				vPos[1] = vertices[vi][1] + noiseCoeff * rnd.nextDouble();
				vPos[2] = vertices[vi][2] + noiseCoeff * rnd.nextDouble();
				v.P = vPos.clone();
				Pn.setToLength(v.P, v.P, 1.0, Pn.EUCLIDEAN);
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
			writeErrorLine(i + "\t" + absErr + "\t" + argErr + "\t" + reErr + "\t" + imErr + "\t" + meshQuality);
		}

	}

}
