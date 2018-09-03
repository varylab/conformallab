package de.varylab.discreteconformal.holomorphicformsexperiments.convergence;

import static de.varylab.discreteconformal.util.LaplaceUtility.calculateCotanWeights;

import java.io.FileWriter;
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
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.adapter.MappedWeightAdapter;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.holomorphicformsexperiments.Utility;
import de.varylab.discreteconformal.logging.LoggingUtility;
import de.varylab.discreteconformal.math.ComplexUtility;
import de.varylab.discreteconformal.plugin.HyperellipticCurvePlugin;
import de.varylab.discreteconformal.unwrapper.SphericalNormalizerPETSc;
import de.varylab.discreteconformal.util.DiscreteRiemannUtility;
import de.varylab.discreteconformal.util.DiscreteRiemannUtility.Result;

public class LawsonG2Scattering {

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
		"title(\"g = 2\")\n\n";
	
	public static String filename = "LawsonG2-Scatter-01";
	
	public static void main(String[] args) throws Exception {
		LoggingUtility.initLogging();
		final AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		Complex[] b_raw = new Complex[] {
			new Complex(-0.5, -0.8660254037844386),
			new Complex(0.5, -0.8660254037844386),
			new Complex(1.0, 0.0),
			new Complex(0.5, 0.8660254037844386),
			new Complex(-0.5, 0.8660254037844386),
			new Complex(-1.0, 0.0)
		};
		final Complex[] b = normalizeBranchPoints(a, b_raw);
		final int count = 1;
		final Random rnd = new Random();
		
		// Clustering Random
		final CoHDS[] CR_rS = new CoHDS[count];
		final ComplexMatrix[] CR_rPo = new ComplexMatrix[count];
		final ComplexMatrix[] CR_rP = new ComplexMatrix[count];
		final double[] CR_rPnorm = new double[count];
		final double[] CR_maxLengths = new double[count];
		final double[] CR_minAngle = new double[count];
		Thread CR_job = new Thread(new Runnable() {
			@Override
			public void run() {
				Formatter f;
				try {
					FileWriter out = new FileWriter(filename + "-CR.dat", true);
					f = new Formatter(out);
				} catch (Exception e) { System.err.println(e); return; }
				while (true) {
					try {
						int numPoints = rnd.nextInt(2000);
						System.out.println("CR " + numPoints);
						calculate(
							count, // number of iterations
							a, // adapters
							rnd, // random number generator
							b, // branch data
							true, // use clustering
							new int[] {0}, // homogeneous points
							new int[] {numPoints}, // clustered points
							false, // use fibonacci for clustering
							0, // equalizer iterations to apply
							CR_rPo, // original period matrix output
							CR_rP, // period matrix output
							CR_rS, // surface output
							CR_rPnorm, // norm of period matrix output
							CR_maxLengths, // effective maximum edge length
							CR_minAngle // minimum triangle angle
						);
						f.format("%1$s, %2$s, %3$s, %4$d\n", CR_rPnorm[0], CR_maxLengths[0], CR_minAngle[0], CR_rS[0].numVertices());
						f.flush();
					} catch (Exception e) {
						System.err.println(e);
					}
				}
			}
		});
		
		// Clustering Fibonacci
		final CoHDS[] CF_rS = new CoHDS[count];
		final ComplexMatrix[] CF_rPo = new ComplexMatrix[count];
		final ComplexMatrix[] CF_rP = new ComplexMatrix[count];
		final double[] CF_rPnorm = new double[count];
		final double[] CF_maxLengths = new double[count];
		final double[] CF_minAngle = new double[count];
		Thread CF_job = new Thread(new Runnable() {
			@Override
			public void run() {
				Formatter f;
				try {
					FileWriter out = new FileWriter(filename + "-CF.dat", true);
					f = new Formatter(out);
				} catch (Exception e) { System.err.println(e); return; }
				while (true) {
					try {
						int numPoints = rnd.nextInt(2000);
						System.out.println("CF " + numPoints);				
						calculate(
							count, // number of iterations
							a, // adapters
							rnd, // random number generator
							b, // branch data
							true, // use clustering
							new int[] {0}, // homogeneous points
							new int[] {numPoints}, // clustered points
							true, // use fibonacci for clustering
							0, // equalizer iterations to apply
							CF_rPo, // original period matrix output
							CF_rP, // period matrix output
							CF_rS, // surface output
							CF_rPnorm, // norm of period matrix output
							CF_maxLengths, // effective maximum edge length
							CF_minAngle // minimum triangle angle
						);
						f.format("%1$s, %2$s, %3$s, %4$d\n", CF_rPnorm[0], CF_maxLengths[0], CF_minAngle[0], CF_rS[0].numVertices());
						f.flush();
					} catch (Exception e) {
						System.err.println(e);
					}
				}
			}
		});
		
		// Homogeneous Random
		final CoHDS[] HR_rS = new CoHDS[count];
		final ComplexMatrix[] HR_rPo = new ComplexMatrix[count];
		final ComplexMatrix[] HR_rP = new ComplexMatrix[count];
		final double[] HR_rPnorm = new double[count];
		final double[] HR_maxLengths = new double[count];
		final double[] HR_minAngle = new double[count];
		Thread HR_job = new Thread(new Runnable() {
			@Override
			public void run() {
				Formatter f;
				try {
					FileWriter out = new FileWriter(filename + "-HR.dat", true);
					f = new Formatter(out);
				} catch (Exception e) { System.err.println(e); return; }
				while (true) {
					try {
						int numPoints = rnd.nextInt(8000);
						System.out.println("HR " + numPoints);					
						calculate(
							count, // number of iterations
							a, // adapters
							rnd, // random number generator
							b, // branch data
							false, // use clustering
							new int[] {numPoints}, // homogeneous points
							new int[] {0}, // clustered points
							false, // use fibonacci for clustering
							0, // equalizer iterations to apply
							HR_rPo, // original period matrix output
							HR_rP, // period matrix output
							HR_rS, // surface output
							HR_rPnorm, // norm of period matrix output
							HR_maxLengths, // effective maximum edge length
							HR_minAngle // minimum triangle angle
						);
						f.format("%1$s, %2$s, %3$s, %4$d\n", HR_rPnorm[0], HR_maxLengths[0], HR_minAngle[0], HR_rS[0].numVertices());
						f.flush();
					} catch (Exception e) {
						System.err.println(e);
					}
				}				
			}
		});
		
		// Homogeneous Fibonacci
		final CoHDS[] HF_rS = new CoHDS[count];
		final ComplexMatrix[] HF_rPo = new ComplexMatrix[count];
		final ComplexMatrix[] HF_rP = new ComplexMatrix[count];
		final double[] HF_rPnorm = new double[count];
		final double[] HF_maxLengths = new double[count];
		final double[] HF_minAngle = new double[count];
		Thread HF_job = new Thread(new Runnable() {
			@Override
			public void run() {
				Formatter f;
				try {
					FileWriter out = new FileWriter(filename + "-HF.dat", true);
					f = new Formatter(out);
				} catch (Exception e) { System.err.println(e); return; }
				while (true) {
					try {
						int numPoints = rnd.nextInt(8000);
						System.out.println("HF " + numPoints);					
						calculate(
							count, // number of iterations
							a, // adapters
							rnd, // random number generator
							b, // branch data
							false, // use clustering
							new int[] {numPoints}, // homogeneous points
							new int[] {0}, // clustered points
							true, // use fibonacci for clustering
							0, // equalizer iterations to apply
							HF_rPo, // original period matrix output
							HF_rP, // period matrix output
							HF_rS, // surface output
							HF_rPnorm, // norm of period matrix output
							HF_maxLengths, // effective maximum edge length
							HF_minAngle // minimum triangle angle
						);	
						f.format("%1$s, %2$s, %3$s, %4$d\n", HF_rPnorm[0], HF_maxLengths[0], HF_minAngle[0], HF_rS[0].numVertices());
						f.flush();
					} catch (Exception e) {
						System.err.println(e);
					}
				}							
			}
		});
		
		// start calculation
		CR_job.start();
		CF_job.start();
		HR_job.start();
		HF_job.start();
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
		ComplexMatrix[] rPo,
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
				rPo[i] = r.periodMatrix_original;
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