package de.varylab.discreteconformal.convergence;

import static de.varylab.discreteconformal.util.DiscreteEllipticUtility.calculateHalfPeriodRatio;
import static de.varylab.discreteconformal.util.DiscreteEllipticUtility.generateEllipticImage;
import static java.lang.Math.signum;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
import de.jreality.math.Pn;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class ConvergenceNoise extends ConvergenceSeries {

	private double
		noiseCoeff = 1E-8;
	private int
		numIterations = 1;
	
	@Override
	protected OptionSet configureAndParseOptions(OptionParser p, String... args) {
		OptionSpec<Integer> numIterationsSpec = p.accepts("IT", "Number of iterations").withRequiredArg().ofType(Integer.class).defaultsTo(1);
		OptionSpec<Double> noiseCoefficientSpec = p.accepts("NC", "Noise coefficient").withRequiredArg().ofType(Double.class).defaultsTo(1E-8);
		OptionSet opts = p.parse(args);
		numIterations = numIterationsSpec.value(opts);
		noiseCoeff = noiseCoefficientSpec.value(opts);
		return opts;
	}
	
	@Override
	protected void perform() throws Exception {
		writeComment("index[1], absErr[2], argErr[3], reErr[4], imErr[5]");
		
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
			CoVertex cutRoot = hds.getVertex(branchIndices[0]);
			CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<>();
			try {
				Set<CoEdge> glueSet = new HashSet<CoEdge>();
				Map<CoVertex, CoVertex> involution = generateEllipticImage(hds, 0, glueSet, branchIndices);
				if (!cutRoot.isValid()) cutRoot = involution.get(cutRoot);
				tau = calculateHalfPeriodRatio(hds, cutRoot, 1E-8, cutInfo);
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
			writeData(i + "\t" + absErr + "\t" + argErr + "\t" + reErr + "\t" + imErr);
		}

	}

}
