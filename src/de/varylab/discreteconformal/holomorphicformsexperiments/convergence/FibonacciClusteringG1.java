package de.varylab.discreteconformal.holomorphicformsexperiments.convergence;

import static de.varylab.discreteconformal.util.LaplaceUtility.calculateCotanWeights;

import java.util.Arrays;
import java.util.Formatter;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import de.jtem.blas.ComplexMatrix;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.generic.UndirectedEdgeIndex;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.generic.Position3d;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.io.HalfedgeIO;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.surface.BranchPoint;
import de.jtem.riemann.theta.SiegelReduction;
import de.varylab.discreteconformal.adapter.MappedWeightAdapter;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.holomorphicformsexperiments.Utility;
import de.varylab.discreteconformal.logging.LoggingUtility;
import de.varylab.discreteconformal.math.ComplexUtility;
import de.varylab.discreteconformal.plugin.HyperellipticCurvePlugin;
import de.varylab.discreteconformal.plugin.hyperelliptic.Curve;
import de.varylab.discreteconformal.unwrapper.SphericalNormalizerPETSc;
import de.varylab.discreteconformal.util.DiscreteRiemannUtility;
import de.varylab.discreteconformal.util.DiscreteRiemannUtility.Result;

public class FibonacciClusteringG1 {

	public static String rTemplate = 
		"ref = %1$s\n" + 
		"\n" + 
		"clusteringRandom <- abs(c(%2$s) - ref)\n" +
		"clusteringRandom_range <- c(%3$s)\n" + 
		"clusteringFibonacci <- abs(c(%4$s) - ref)\n" +
		"clusteringFibonacci_range <- c(%5$s)\n" +
		"homogeneousRandom <- abs(c(%6$s) - ref)\n" + 
		"homogeneousRandom_range <- c(%7$s)\n" + 
		"homogeneousFibonacci <- abs(c(%8$s) - ref)\n" +
		"homogeneousFibonacci_range <- c(%9$s)\n" +
		"\n" +
		"xlim = range(clusteringRandom_range, clusteringFibonacci_range, homogeneousRandom_range, homogeneousFibonacci_range)\n" + 
		"ylim = range(clusteringRandom, clusteringFibonacci, homogeneousRandom, homogeneousFibonacci)\n" + 
		"\n" + 
		"plot(\n" + 
		"	y=clusteringRandom, ylab=\"Error\", ylim=ylim,\n" + 
		"	x=clusteringRandom_range, xlab=\"Max Edge Length\", xlim=xlim,\n" + 
		"	log=\"xy\", type=\"o\", pch=22, lty=1, col=1\n" + 
		")\n" + 
		"lines(y=clusteringFibonacci, x=clusteringFibonacci_range, type=\"o\", pch=22, lty=1, col=2)\n" + 
		"lines(y=homogeneousRandom, x=homogeneousRandom_range, type=\"o\", pch=22, lty=1, col=3)\n" + 
		"lines(y=homogeneousFibonacci, x=homogeneousFibonacci_range, type=\"o\", pch=22, lty=1, col=4)\n" + 
		"\n" + 
		"legend(\"topright\", c(\"Clustering Random\",\"Clustering Fibonacci\",\"Homogeneous Random\",\"Homogeneous Fibonacci\"), col=c(1,2,3,4), lty=1);\n" + 
		"title(\"General Torus\")";
	public static String filename = "FibonacciClusteringG1-01";
	
	public static void main(String[] args) throws Exception {
		LoggingUtility.initLogging();
		final AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		Complex[] b_raw = new Complex[] {
			new Complex(0.5, 0.4),
			new Complex(-0.3, 0.2),
			new Complex(-0.1, -0.0),
			new Complex(0.1, -0.2)
		};
		final Complex[] b = normalizeBranchPoints(a, b_raw);
		
		
		BranchPoint[] bb = new BranchPoint[b.length];
		for (int i = 0; i < b.length; i++) bb[i] = new BranchPoint(b[i]);
		Curve C = new Curve(bb);
		C.setEps(1E-15);
		ComplexMatrix P = C.getPeriodMatrix();
		SiegelReduction siegel = new SiegelReduction(P);
		P = siegel.getReducedPeriodMatrix();

		final int count = 4;
		final Random rnd = new Random(2);
		
		// Clustering Random
		final CoHDS[] CR_rS = new CoHDS[count];
		final ComplexMatrix[] CR_rP = new ComplexMatrix[count];
		final double[] CR_rPnorm = new double[count];
		final double[] CR_maxLengths = new double[count];
		final double[] CR_minAngle = new double[count];
		Thread CR_job = new Thread(new Runnable() {
			@Override
			public void run() {
				FibonacciClusteringG1.calculate(
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
			}
		});
		
		// Clustering Fibonacci
		final CoHDS[] CF_rS = new CoHDS[count];
		final ComplexMatrix[] CF_rP = new ComplexMatrix[count];
		final double[] CF_rPnorm = new double[count];
		final double[] CF_maxLengths = new double[count];
		final double[] CF_minAngle = new double[count];
		Thread CF_job = new Thread(new Runnable() {
			@Override
			public void run() {
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
			}
		});
		
		// Homogeneous Random
		final CoHDS[] HR_rS = new CoHDS[count];
		final ComplexMatrix[] HR_rP = new ComplexMatrix[count];
		final double[] HR_rPnorm = new double[count];
		final double[] HR_maxLengths = new double[count];
		final double[] HR_minAngle = new double[count];
		Thread HR_job = new Thread(new Runnable() {
			@Override
			public void run() {
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
			}
		});
		
		// Homogeneous Fibonacci
		final CoHDS[] HF_rS = new CoHDS[count];
		final ComplexMatrix[] HF_rP = new ComplexMatrix[count];
		final double[] HF_rPnorm = new double[count];
		final double[] HF_maxLengths = new double[count];
		final double[] HF_minAngle = new double[count];
		Thread HF_job = new Thread(new Runnable() {
			@Override
			public void run() {
				calculate(
					count, // number of iterations
					a, // adapters
					rnd, // random number generator
					b, // branch data
					false, // use clustering
					new int[] {40, 400, 4000, 8000, 16000, 32000}, // homogeneous points
					new int[] {0, 0, 0, 0, 0, 0}, // clustered points
					true, // use fibonacci for clustering
					0, // equalizer iterations to apply
					HF_rP, // period matrix output
					HF_rS, // surface output
					HF_rPnorm, // norm of period matrix output
					HF_maxLengths, // effective maximum edge length
					HF_minAngle // minimum triangle angle
				);		
			}
		});
		
		// start calculation
		CR_job.start();
		CF_job.start();
		HR_job.start();
		HF_job.start();
		
		// wait for results
		CR_job.join();
		CF_job.join();
		HR_job.join();
		HF_job.join();
		
		Formatter out = new Formatter(FibonacciClusteringG1.filename + ".r");
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
		out.format("# Branch Data: %1$s\n", Arrays.toString(b_raw));
		out.format("#  REF  : %1$s\n", P);
		out.format("# |REF| : %1$s\n", P.normSqr());
		out.format("# Clustering Random ---------\n");
		for (int i = 0; i < count; i++) {
			out.format("#  S[%2$d] : %1$s\n", CR_rS[i], i);
			out.format("#  P[%2$d] : %1$s\n", CR_rP[i], i);
			out.format("# |P[%2$d]|: %1$s\n", CR_rP[i].normSqr(), i);
			out.format("# |L[%2$d]|: %1$s\n", CR_maxLengths[i], i);
			out.format("# |a[%2$d]|: %1$s\n", CR_minAngle[i], i);
			HalfedgeIO.writeOBJ(CR_rS[i], a, FibonacciClusteringG1.filename + "_CR_" + i + ".obj");
		}
		out.format("# Clustering Fibonacci ------\n");
		for (int i = 0; i < count; i++) {
			out.format("#  S[%2$d] : %1$s\n", CF_rS[i], i);
			out.format("#  P[%2$d] : %1$s\n", CF_rP[i], i);
			out.format("# |P[%2$d]|: %1$s\n", CF_rP[i].normSqr(), i);
			out.format("# |L[%2$d]|: %1$s\n", CF_maxLengths[i], i);
			out.format("# |a[%2$d]|: %1$s\n", CF_minAngle[i], i);
			HalfedgeIO.writeOBJ(CF_rS[i], a, FibonacciClusteringG1.filename + "_CF_" + i + ".obj");
		}
		out.format("# Homogeneous Random --------\n");
		for (int i = 0; i < count; i++) {
			out.format("#  S[%2$d] : %1$s\n", HR_rS[i], i);
			out.format("#  P[%2$d] : %1$s\n", HR_rP[i], i);
			out.format("# |P[%2$d]|: %1$s\n", HR_rP[i].normSqr(), i);
			out.format("# |L[%2$d]|: %1$s\n", HR_maxLengths[i], i);
			out.format("# |a[%2$d]|: %1$s\n", HR_minAngle[i], i);
			HalfedgeIO.writeOBJ(HR_rS[i], a, FibonacciClusteringG1.filename + "_HR_" + i + ".obj");
		}
		out.format("# Homogeneous Fibonacci -----\n");
		for (int i = 0; i < count; i++) {
			out.format("#  S[%2$d] : %1$s\n", HF_rS[i], i);
			out.format("#  P[%2$d] : %1$s\n", HF_rP[i], i);
			out.format("# |P[%2$d]|: %1$s\n", HF_rP[i].normSqr(), i);
			out.format("# |L[%2$d]|: %1$s\n", HF_maxLengths[i], i);
			out.format("# |a[%2$d]|: %1$s\n", HF_minAngle[i], i);
			HalfedgeIO.writeOBJ(HF_rS[i], a, FibonacciClusteringG1.filename + "_HF_" + i + ".obj");
		}
		out.close();
	}

	private static Complex[] normalizeBranchPoints(final AdapterSet a, Complex[] b_raw) {
		CoHDS branchHDS = new CoHDS();
		for (Complex bp : b_raw) {
			double[] pos = ComplexUtility.inverseStereographic(bp);
			CoVertex bv = branchHDS.addNewVertex();
			a.set(Position.class, bv, pos);
		}
		SphericalNormalizerPETSc.normalize(branchHDS, a, Position4d.class, Position.class);
		final Complex[] b = new Complex[b_raw.length];
		for (CoVertex v : branchHDS.getVertices()) {
			double[] pos = a.getD(Position3d.class, v);
			b[v.getIndex()] = ComplexUtility.stereographic(pos);
		}
		return b;
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

	public static void calculate(
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
			try {
				Set<CoVertex> branchVertices = new HashSet<>();
				CoHDS S = HyperellipticCurvePlugin.generateCurve(
					b, false, numextra[i], numextrabranch[i], useFibonacciPoints, numEqualizerIterations, rnd, a, branchVertices
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
			} catch (Exception e) {
				System.err.println(e);
			}
		}
	}
	
}
