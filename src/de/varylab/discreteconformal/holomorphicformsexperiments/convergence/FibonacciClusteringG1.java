package de.varylab.discreteconformal.holomorphicformsexperiments.convergence;

import static de.varylab.discreteconformal.util.LaplaceUtility.calculateCotanWeights;

import java.util.Formatter;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import de.jtem.blas.ComplexMatrix;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.generic.UndirectedEdgeIndex;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.surface.BranchPoint;
import de.jtem.riemann.theta.SiegelReduction;
import de.varylab.discreteconformal.adapter.MappedWeightAdapter;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.holomorphicformsexperiments.Utility;
import de.varylab.discreteconformal.logging.LoggingUtility;
import de.varylab.discreteconformal.plugin.HyperellipticCurvePlugin;
import de.varylab.discreteconformal.plugin.hyperelliptic.Curve;
import de.varylab.discreteconformal.util.DiscreteRiemannUtility;
import de.varylab.discreteconformal.util.DiscreteRiemannUtility.Result;

public class FibonacciClusteringG1 {

	public static String rTemplate = 
		"ref = %1$s\n" + 
		"\n" + 
		"clusteringRandom <- abs(c(%2$s) - ref)\n" +
		"clusteringRandom_range <- c(%3$s)\n" + 
		"clusteringFibonacci <- abs(c(%4$s) - ref)\n" +
		"clusteringFibonacci_range <- abs(c(%5$s) - ref)\n" +
		"homogeneousRandom <- abs(c(%6$s) - ref)\n" + 
		"homogeneousRandom_range <- abs(c(%7$s) - ref)\n" + 
		"homogeneousFibonacci <- abs(c(%8$s) - ref)\n" +
		"homogeneousFibonacci_range <- abs(c(%9$s) - ref)\n" +
		"\n" + 
		"plot(y=clusteringRandom, ylab=\"Error\", x=clusteringRandom_range, xlab=\"Max Edge Length\", log=\"xy\", type=\"o\", pch=22, lty=1, col=1\n" +
		"lines(y=clusteringFibonacci, x=clusteringFibonacci_range, type=\"o\", pch=22, lty=1, col=2, ylim=g_range)\n" +
		"lines(y=homogeneousRandom, x=homogeneousRandom_range, type=\"o\", pch=22, lty=1, col=3, ylim=g_range)\n" + 
		"lines(y=homogeneousFibonacci, x=homogeneousFibonacci_range, type=\"o\", pch=22, lty=1, col=4, ylim=g_range)\n" + 
		"\n" + 
		"legend(\"topright\", c(\"Clustering Random\",\"Clustering Fibonacci\",\"Homogeneous Random\",\"Homogeneous Fibonacci\"), col=c(1,2,3,4), lty=1);\n" + 
		"title(\"General Torus\")\n" +
		"\n";
	
	public static void main(String[] args) throws Exception {
		LoggingUtility.initLogging();
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		Random rnd = new Random(4);
		Complex[] b = new Complex[] {
			new Complex(0.5, 0.4),
			new Complex(-0.3, 0.2),
			new Complex(-0.1, -0.0),
			new Complex(0.1, -0.2)
		};
		BranchPoint[] bb = new BranchPoint[b.length];
		for (int i = 0; i < b.length; i++) bb[i] = new BranchPoint(b[i]);
		Curve C = new Curve(bb);
		C.setEps(1E-15);
		ComplexMatrix P = C.getPeriodMatrix();
		SiegelReduction siegel = new SiegelReduction(P);
		P = siegel.getReducedPeriodMatrix();

		int count = 2;
		
		// Clustering Random
		CoHDS[] CR_rS = new CoHDS[count];
		ComplexMatrix[] CR_rP = new ComplexMatrix[count];
		double[] CR_rPnorm = new double[count];
		double[] CR_maxLengths = new double[count];
		double[] CR_minAngle = new double[count];
		calculate(
			count, // number of iterations
			a, // adapters
			rnd, // random number generator
			b, // branch data
			true, // use clustering
			new int[] {0, 0, 0, 0, 0, 0}, // homogeneous points
			new int[] {10, 100, 1000, 2000, 4000, 8000}, // clustered points
			false, // use fibonacci for clustering
			0, // equalizer iterations to apply
			CR_rP, // period matrix output
			CR_rS, // surface output
			CR_rPnorm, // norm of period matrix output
			CR_maxLengths, // effective maximum edge length
			CR_minAngle // minimum triangle angle
		);
		
		// Clustering Fibonacci
		CoHDS[] CF_rS = new CoHDS[count];
		ComplexMatrix[] CF_rP = new ComplexMatrix[count];
		double[] CF_rPnorm = new double[count];
		double[] CF_maxLengths = new double[count];
		double[] CF_minAngle = new double[count];
		calculate(
			count, // number of iterations
			a, // adapters
			rnd, // random number generator
			b, // branch data
			true, // use clustering
			new int[] {0, 0, 0, 0, 0, 0}, // homogeneous points
			new int[] {10, 100, 1000, 2000, 4000, 8000}, // clustered points
			true, // use fibonacci for clustering
			0, // equalizer iterations to apply
			CF_rP, // period matrix output
			CF_rS, // surface output
			CF_rPnorm, // norm of period matrix output
			CF_maxLengths, // effective maximum edge length
			CF_minAngle // minimum triangle angle
		);
		
		// Homogeneous Random
		CoHDS[] HR_rS = new CoHDS[count];
		ComplexMatrix[] HR_rP = new ComplexMatrix[count];
		double[] HR_rPnorm = new double[count];
		double[] HR_maxLengths = new double[count];
		double[] HR_minAngle = new double[count];
		calculate(
			count, // number of iterations
			a, // adapters
			rnd, // random number generator
			b, // branch data
			false, // use clustering
			new int[] {40, 400, 4000, 8000, 16000, 32000}, // homogeneous points
			new int[] {0, 0, 0, 0, 0, 0}, // clustered points
			false, // use fibonacci for clustering
			0, // equalizer iterations to apply
			HR_rP, // period matrix output
			HR_rS, // surface output
			HR_rPnorm, // norm of period matrix output
			HR_maxLengths, // effective maximum edge length
			HR_minAngle // minimum triangle angle
		);
		
		// Homogeneous Fibonacci
		CoHDS[] HF_rS = new CoHDS[count];
		ComplexMatrix[] HF_rP = new ComplexMatrix[count];
		double[] HF_rPnorm = new double[count];
		double[] HF_maxLengths = new double[count];
		double[] HF_minAngle = new double[count];
		calculate(
			count, // number of iterations
			a, // adapters
			rnd, // random number generator
			b, // branch data
			false, // use clustering
			new int[] {40, 400, 4000, 8000, 16000, 32000}, // homogeneous points
			new int[] {0, 0, 0, 0, 0, 0}, // clustered points
			false, // use fibonacci for clustering
			0, // equalizer iterations to apply
			HF_rP, // period matrix output
			HF_rS, // surface output
			HF_rPnorm, // norm of period matrix output
			HF_maxLengths, // effective maximum edge length
			HF_minAngle // minimum triangle angle
		);		
		
		String packagePath = FibonacciClusteringG1.class.getName().replace('.', '/');
		String filePath = String.format("src/%1$s.r", packagePath);
		Formatter out = new Formatter(filePath);
		out.format(FibonacciClusteringG1.rTemplate,
			P.normSqr(),
			join(CR_rPnorm),
			join(CR_maxLengths),
			join(CF_rPnorm),
			join(CF_maxLengths),
			join(HR_rPnorm),
			join(HR_maxLengths),
			join(HF_rPnorm),
			join(HF_maxLengths)
		);
		
		out.format("# ---------------------------\n");
		out.format("#  REF  : %1$s\n", P);
		out.format("# |REF| : %1$s\n", P.normSqr());
		out.format("# Clustering Random ---------\n");
		for (int i = 0; i < count; i++) {
			out.format("#  S[%2$d] : %1$s\n", CR_rS[i], i);
			out.format("#  P[%2$d] : %1$s\n", CR_rP[i], i);
			out.format("# |P[%2$d]|: %1$s\n", CR_rP[i].normSqr(), i);
			out.format("# |L[%2$d]|: %1$s\n", CR_maxLengths[i], i);
			out.format("# |a[%2$d]|: %1$s\n", CR_minAngle[i], i);
		}
		out.format("# Clustering Fibonacci ------\n");
		for (int i = 0; i < count; i++) {
			out.format("#  S[%2$d] : %1$s\n", CF_rS[i], i);
			out.format("#  P[%2$d] : %1$s\n", CF_rP[i], i);
			out.format("# |P[%2$d]|: %1$s\n", CF_rP[i].normSqr(), i);
			out.format("# |L[%2$d]|: %1$s\n", CF_maxLengths[i], i);
			out.format("# |a[%2$d]|: %1$s\n", CF_minAngle[i], i);
		}
		out.format("# Homogeneous Random --------\n");
		for (int i = 0; i < count; i++) {
			out.format("#  S[%2$d] : %1$s\n", HR_rS[i], i);
			out.format("#  P[%2$d] : %1$s\n", HR_rP[i], i);
			out.format("# |P[%2$d]|: %1$s\n", HR_rP[i].normSqr(), i);
			out.format("# |L[%2$d]|: %1$s\n", HR_maxLengths[i], i);
			out.format("# |a[%2$d]|: %1$s\n", HR_minAngle[i], i);
		}
		out.format("# Homogeneous Fibonacci -----\n");
		for (int i = 0; i < count; i++) {
			out.format("#  S[%2$d] : %1$s\n", HF_rS[i], i);
			out.format("#  P[%2$d] : %1$s\n", HF_rP[i], i);
			out.format("# |P[%2$d]|: %1$s\n", HF_rP[i].normSqr(), i);
			out.format("# |L[%2$d]|: %1$s\n", HF_maxLengths[i], i);
			out.format("# |a[%2$d]|: %1$s\n", HF_minAngle[i], i);
		}
		out.close();
	}

	private static String join(double[] CR_rPnorm) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < CR_rPnorm.length; i++) {
			double d = CR_rPnorm[i];
			if (i != 0) { sb.append(','); }
			sb.append(d);
		}
		return sb.toString();
	}

	private static void calculate(
		int count,
		AdapterSet a,
		Random rnd,
		Complex[] b,
		boolean clustering,
		int[] numextra,
		int[] numextrabranch,
		boolean useFibonacciPoints,
		int numEqualizerIterations,
		ComplexMatrix[] rP,
		CoHDS[] rS,
		double[] rPnorm,
		double[] maxLengths,
		double[] minAngle
	) {
		for (int i = 0; i < count; i++) {
			Set<CoVertex> branchVertices = new HashSet<>();
			CoHDS S = HyperellipticCurvePlugin.generateCurve(
				b, numextra[i], numextrabranch[i], useFibonacciPoints, numEqualizerIterations, rnd, a, branchVertices
			);
			MappedWeightAdapter cotanWeights = calculateCotanWeights(S, a);
			AdapterSet aa = new AdapterSet(a);
			aa.add(new UndirectedEdgeIndex());
			aa.add(cotanWeights);
			Result r = DiscreteRiemannUtility.getHolomorphicFormsAndPeriodMatrix(S, aa, null);
			rP[i] = r.periodMatrix;
			rS[i] = S;
			rPnorm[i] = r.periodMatrix.normSqr();
			maxLengths[i] = Utility.calculateLargestEdgeLength(S, a, branchVertices, clustering);
			minAngle[i] = Utility.calculateSmallestAngle(S, a);
		}
	}
	
}
