package de.varylab.discreteconformal.util;

import java.util.List;

import cern.colt.matrix.Norm;
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.solver.DefaultDoubleIterationMonitor;
import cern.colt.matrix.tdouble.algo.solver.DoubleBiCGstab;
import cern.colt.matrix.tdouble.algo.solver.DoubleGMRES;
import cern.colt.matrix.tdouble.algo.solver.DoubleIterationReporter;
import cern.colt.matrix.tdouble.algo.solver.DoubleIterativeSolver;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.algorithm.triangulation.MappedLengthAdapter;
import de.varylab.discreteconformal.util.Search.WeightAdapter;


/**
 * Class to calculate holomorphic differentials in the sense of mercat.
 * 
 * By convention all the forms are saved as rows and all the cycles as columns.
 * (For private methods we always use colt matrices.)
 * 
 * @author knoeppel
 * 
 */
public class DiscreteHolomorphicFormUtility {

	private static DoubleIterationReporter reporter = new ColtIterationReporterImpl();
	
	private static DenseDoubleAlgebra dalgebra = new DenseDoubleAlgebra();
	
	/**
	 * Returns a basis of holomorphic differentials on the surface, which is
	 * normalized with respect to the given homology basis (canonical).
	 * The forms are returned as a two-array of matrices, the first entry of
	 * which belongs to the real and the second to the imaginary part. The
	 * surface has to be Delaunay triangulated. The forms are stored as rows.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param canonicalHomologyBasis
	 * @param adapters
	 * @param la
	 * @param wa
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix2D[] getHolomorphicFormsOnPrimalMesh(
				HalfEdgeDataStructure<V, E, F> delaunay,
				List<List<E>> canonicalHomologyBasis, AdapterSet adapters,
				MappedLengthAdapter la, WeightAdapter<E> wa) {

		// the forms are defined on the positive oriented edges
		int numPosEdges = delaunay.numEdges() / 2;
		
		// g is simply the genus of the surface
		int g = canonicalHomologyBasis.size() / 2;

		// Get the harmonic differentials on the surface and its dual. The
		// format of the matrices is 2g*numEdges
		DoubleMatrix2D dh = DiscreteHarmonicFormUtility
				.getHarmonicFormsOfPrimalMesh(delaunay, canonicalHomologyBasis,
						adapters, la);
		DoubleMatrix2D dhStar = DualityUtility.getDualOfPrimalForms(delaunay,
				adapters, dh);

		// to normalize the differentials we need the a-periods and its dual
		// cycles
		List<List<E>> acycles = CanonicalBasisUtility
				.getACycles(canonicalHomologyBasis);
		List<List<E>> dualacycles = DualityUtility.getDualPaths(delaunay,
				acycles);

		// write cycles to matrices
		DoubleMatrix2D A = EdgeUtility.cyclesToMatrix(adapters, delaunay,
				acycles);
		DoubleMatrix2D dualA = EdgeUtility.cyclesToMatrix(adapters, delaunay,
				dualacycles);

		// an array to save the forms: The real part corresponds to the
		// differential on the primal mesh and the imaginary part to the
		// differentials on the dual mesh.
		DoubleMatrix2D[] OMEGA = new DoubleMatrix2D[2];
		for (int i = 0; i < OMEGA.length; i++)
			OMEGA[i] = DoubleFactory2D.dense.make(g, numPosEdges);

		// The holomorphic forms omega are linear combinations of 2g harmonic
		// forms dh_j plus i times its duals (star)dh_j: omega= sum(x_j*(dh_j+
		// i*(star)dh_j)). They are normalized iff they satisfy a linear system.
		// A is the coefficient matrix.
		DoubleMatrix2D M = DiscreteHarmonicFormUtility
				.getHarmonicPeriodsMatrix(dh, dhStar, A, dualA);

		// print(M, 2);

		// The number of harmonic forms on the surface is 2g, so we need 2g
		// coefficients.
		DoubleMatrix1D x = DoubleFactory1D.dense.make(2 * g);
		
		// For each a-cycle find the holomorphic form which is 1 along this
		// cycle, i.e. gives one for the primal cycle and 0 for its dual cycle. 
		for (int i = 0; i < g; i++) {
			
			// set up the conditions 
			DoubleMatrix1D bc = DoubleFactory1D.dense.make(2 * g);
			bc.set(g+i, 2*Math.PI);
			
			// solve the system
			solve(M, x, bc);

			// Build linear combinations of the harmonic differentials
			// using the obtained coefficients. 
			DoubleMatrix1D omega= dalgebra.mult(dalgebra.transpose(dh), x);
			DoubleMatrix1D omegaStar= dalgebra.mult(dalgebra.transpose(dhStar), x);
			
			for (int j = 0; j < numPosEdges; j++) {
				OMEGA[0].set(i, j, omega.get(j));
				OMEGA[1].set(i, j, omegaStar.get(j));
			}
		}
		
		adapters.remove(la);

		return OMEGA;
	}
	
	/**
	 * Returns a basis of holomorphic differentials on the dual of the surface,
	 * which is normalized with respect to the given homology basis (canonical).
	 * The forms are returned as a two-array of matrices, the first entry of
	 * which belongs to the real and the second to the imaginary part. The
	 * surface has to be Delaunay triangulated. The forms are stored as rows.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param canonicalHomologyBasis
	 * @param adapters
	 * @param la
	 * @param wa
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix2D[] getHolomorphicFormsOnDualMesh(
				HalfEdgeDataStructure<V, E, F> delaunay,
				List<List<E>> canonicalHomologyBasis, AdapterSet adapters,
				MappedLengthAdapter la, WeightAdapter<E> wa) {

		// the forms are defined on the positive oriented edges
		int numPosEdges = delaunay.numEdges() / 2;

		// g is simply the genus of the surface
		int g = canonicalHomologyBasis.size() / 2;

		List<List<E>> dualBasis = DualityUtility.getDualPaths(delaunay,
				canonicalHomologyBasis);

		// Get the harmonic differentials on the surface and its dual. The
		// format of the matrices is 2g*numEdges
		DoubleMatrix2D dhStar = DiscreteHarmonicFormUtility
				.getHarmonicFormsOfDualMesh(delaunay, dualBasis,
						adapters, la);
		DoubleMatrix2D dhStarStar = DualityUtility.getDualOfDualForms(delaunay,
				adapters, dhStar);

		// to normalize the differentials we need the a-periods and its dual
		// cycles
		List<List<E>> dualacycles = CanonicalBasisUtility.getACycles(dualBasis);
		List<List<E>> primalacycles = CanonicalBasisUtility
				.getACycles(canonicalHomologyBasis);

		// write cycles to matrices
		DoubleMatrix2D dualA = EdgeUtility.cyclesToMatrix(adapters, delaunay,
				dualacycles);
		DoubleMatrix2D primalA = EdgeUtility.cyclesToMatrix(adapters, delaunay,
				primalacycles);

		// an array to save the forms: The real part corresponds to the
		// differential on the primal mesh and the imaginary part to the
		// differentials on the dual mesh.
		DoubleMatrix2D[] OMEGA = new DoubleMatrix2D[2];
		for (int i = 0; i < OMEGA.length; i++)
			OMEGA[i] = DoubleFactory2D.dense.make(g, numPosEdges);

		// The holomorphic forms omega are linear combinations of 2g harmonic
		// forms dh_j plus i times its duals (star)dh_j: omega= sum(x_j*(dh_j+
		// i*(star)dh_j)). They are normalized iff they satisfy a linear system.
		// A is the coefficient matrix.
		DoubleMatrix2D M = DiscreteHarmonicFormUtility
				.getHarmonicPeriodsMatrix(dhStar, dhStarStar, dualA, primalA);

		// The number of harmonic forms on the surface is 2g, so we need 2g
		// coefficients.
		DoubleMatrix1D x = DoubleFactory1D.dense.make(2 * g);

		// For each a-cycle find the holomorphic form which is 1 along this
		// cycle, i.e. gives one for the primal cycle and 0 for its dual cycle.
		for (int i = 0; i < g; i++) {

			// set up the conditions
			DoubleMatrix1D bc = DoubleFactory1D.dense.make(2 * g);
			bc.set(g+i, 2 * Math.PI);

			// solve the system
			solve(M, x, bc);

			// Build linear combinations of the harmonic differentials
			// using the obtained coefficients.
			DoubleMatrix1D omega = dalgebra.mult(dalgebra.transpose(dhStar), x);
			DoubleMatrix1D omegaStar = dalgebra.mult(
					dalgebra.transpose(dhStarStar), x);

			for (int j = 0; j < numPosEdges; j++) {
				OMEGA[0].set(i, j, omega.get(j));
				OMEGA[1].set(i, j, omegaStar.get(j));
			}
		}

		adapters.remove(la);

		return OMEGA;
	}
	
	private static double eps = 1E-20;
	private static int maxIterations= 100000000;
	
	/**
	 * Solves Ax=b and writes the result in the vector x.
	 * 
	 * @param A
	 * @param x
	 * @param b
	 */
	private static void solve(DoubleMatrix2D A, DoubleMatrix1D x,
			DoubleMatrix1D b) {
		
		DoubleIterativeSolver solver;
//		solver = new DoubleGMRES(x);
		solver = new DoubleBiCGstab(x);

		DefaultDoubleIterationMonitor monitor = new DefaultDoubleIterationMonitor();

		// configure monitor
		monitor.setMaxIterations(maxIterations);
		monitor.setAbsoluteTolerance(eps);
		monitor.setRelativeTolerance(eps);
//		monitor.setDivergenceTolerance(1);
		monitor.setNormType(Norm.Infinity);
		monitor.setIterationReporter(reporter);

		solver.setIterationMonitor(monitor);

		try {
			solver.solve(A, b, x);
		} catch (IterativeSolverDoubleNotConvergedException e) {
			System.err
					.println("Iterative solver failed to converge: Couldn't get harmonic function.");
			e.printStackTrace();
		}
		System.err.println();
	}

}
