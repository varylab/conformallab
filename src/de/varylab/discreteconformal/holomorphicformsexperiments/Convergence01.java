package de.varylab.discreteconformal.holomorphicformsexperiments;

import static de.varylab.discreteconformal.util.LaplaceUtility.calculateCotanWeights;

import java.util.Random;

import de.jtem.blas.ComplexMatrix;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.generic.UndirectedEdgeIndex;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.surface.BranchPoint;
import de.jtem.riemann.theta.SiegelReduction;
import de.varylab.discreteconformal.adapter.MappedWeightAdapter;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.logging.LoggingUtility;
import de.varylab.discreteconformal.plugin.HyperellipticCurvePlugin;
import de.varylab.discreteconformal.plugin.hyperelliptic.Curve;
import de.varylab.discreteconformal.util.DiscreteRiemannUtility;
import de.varylab.discreteconformal.util.DiscreteRiemannUtility.Result;

public class Convergence01 {

	public static void main(String[] args) {
		LoggingUtility.initLogging();
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		Random rnd = new Random(1);
		Complex[] b = new Complex[] {
			new Complex(0.5, 0.4),
			new Complex(-0.3, 0.2),
			new Complex(-0.1, -0.0),
			new Complex(0.1, -0.2),
			new Complex(0.3, 0.4),
			new Complex(0.5, 0.6)
		};
		BranchPoint[] bb = new BranchPoint[b.length];
		for (int i = 0; i < b.length; i++) bb[i] = new BranchPoint(b[i]);
		Curve C = new Curve(bb);
		C.setEps(1E-15);
		ComplexMatrix P = C.getPeriodMatrix();
		SiegelReduction siegel = new SiegelReduction(P);
		P = siegel.getReducedPeriodMatrix();
		
		int[] numextra = new int[] {0, 0, 0, 0, 0};
		int[] numextrabranch = new int[] {10, 100, 1000, 2000, 4000};
		int numEqualizerIterations = 5;
		ComplexMatrix[] rP = new ComplexMatrix[numextra.length];
		CoHDS[] rS = new CoHDS[numextra.length];
		for (int i = 0; i < numextra.length; i++) {
			CoHDS S = HyperellipticCurvePlugin.generateCurve(b, numextra[i], numextrabranch[i], numEqualizerIterations, rnd, a, null);
			MappedWeightAdapter cotanWeights = calculateCotanWeights(S, a);
			AdapterSet aa = new AdapterSet(a);
			aa.add(new UndirectedEdgeIndex());
			aa.add(cotanWeights);
			Result r = DiscreteRiemannUtility.getHolomorphicFormsAndPeriodMatrix(S, aa, null);
			rP[i] = r.periodMatrix;
			rS[i] = S;
		}

		System.out.printf(" REF  : %1$s\n", P);
		System.out.printf("|REF| : %1$s\n", P.normSqr());
		for (int i = 0; i < numextra.length; i++) {
			System.out.printf(" S[%2$d] : %1$s\n", rS[i], i);
			System.out.printf(" P[%2$d] : %1$s\n", rP[i], i);
			System.out.printf("|P[%2$d]|: %1$s\n", rP[i].normSqr(), i);
		}
	}
	
}
