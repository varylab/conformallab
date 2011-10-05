package de.varylab.discreteconformal.util;

import static de.varylab.discreteconformal.util.CanonicalBasisUtility.getCanonicalHomologyBasis;
import static de.varylab.discreteconformal.util.DualityUtility.getDualPaths;

import java.util.List;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.decomposition.SparseDoubleLUDecomposition;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.EdgeIndex;
import de.jtem.halfedgetools.adapter.type.Weight;
import de.jtem.halfedgetools.util.HalfEdgeUtilsExtra;
import de.varylab.discreteconformal.adapter.CotanWeightAdapter;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

/**
 * Class to calculate harmonic differentials corresponding to a cycle of the
 * homology basis.
 * 
 * By convention all the forms are saved as rows and all the cycles as columns.
 * (For private methods we always use colt matrices.)
 * 
 * @author knoeppel
 * 
 */
public class DiscreteHarmonicFormUtility {

//	private static DoubleIterationReporter reporter = new ColtIterationReporterImpl();
	private static DenseDoubleAlgebra dalgebra = new DenseDoubleAlgebra();

	/**
	 * Calculates 2*g harmonic differentials on hds. The weight Adapter is used
	 * to find a basis of the homology that are short with respect to the weight
	 * given by wa. The forms are stored in the rows.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param wa
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double[][] getHarmonicFormsOnPrimalMesh(
		HalfEdgeDataStructure<V, E, F> hds, 
		AdapterSet adapters,
		WeightAdapter<E> wa
	) {
		// Get the homology basis of the surface.
		V rootV = hds.getVertex(0);
		List<List<E>> basis = getCanonicalHomologyBasis(rootV, adapters, wa);

		DoubleMatrix2D dh = getHarmonicFormsOfPrimalMesh(hds, basis, adapters);
		SimpleMatrixPrintUtility.print(
			dalgebra.mult(dh, EdgeUtility.cyclesToMatrix(adapters, hds, basis)), 4
		);

		// use the private method
		return dh.toArray();
	}

	/**
	 * Calculates 2*g harmonic differentials on the dual of hds. The weight
	 * Adapter is used to find a basis of the homology that are short with
	 * respect to the weight given by wa. The forms are stored in the rows.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param wa
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double[][] getHarmonicFormsOnDualMesh(
		HalfEdgeDataStructure<V, E, F> hds, 
		AdapterSet adapters,
		WeightAdapter<E> wa
	) {
		// Get the homology basis of the surface.
		V rootV = hds.getVertex(0);
		List<List<E>> dualbasis = getDualPaths(hds, getCanonicalHomologyBasis(rootV, adapters, wa));
		DoubleMatrix2D dh = getHarmonicFormsOfDualMesh(hds, dualbasis, adapters);

		SimpleMatrixPrintUtility.print(
			dalgebra.mult(dh, EdgeUtility.cyclesToMatrix(adapters, hds, dualbasis)), 4
		);

		// use the private method
		return dh.toArray();
	}

	/**
	 * Returns a matrix A needed to build the normalized holomorphic
	 * differentials. The matrix has format 2g*2g. There are 2g harmonic
	 * differentials dH_j= dh1_j+i*dh2_j. The first g rows are filled with its
	 * a1-periods, the second g rows with its values on the a2-cycles. That
	 * means in the upper g rows of Ax stands the real part and in the lower the
	 * imaginary part of the a-periods of the holomorphic differential
	 * omega=sum(x_j*dH_j).
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param DH1
	 * @param DH2
	 * @param A1
	 * @param A2
	 * @param adapters
	 * @return
	 */
	public static DoubleMatrix2D getHarmonicPeriodsMatrix(DoubleMatrix2D DH1,
			DoubleMatrix2D DH2, DoubleMatrix2D A1, DoubleMatrix2D A2) {

		// The genus equals the number of a-cycles equals the number of columns
		// of a
		int g = A1.columns();

		// get a-period-matrix of the harmonic differentials DH
		DoubleMatrix2D aPeriods1 = dalgebra.mult(DH1, A1);
		DoubleMatrix2D aPeriods2 = dalgebra.mult(DH2, A2);

		// The matrix has to have the format 2g*2g.
		DoubleMatrix2D M = DoubleFactory2D.sparse.make(2 * g, 2 * g);

		// for each cycle and each differential
		for (int i = 0; i < g; i++) {
			for (int j = 0; j < 2 * g; j++) {
				// fill the upper g rows with a-periods of the primal mesh
				M.setQuick(i, j, aPeriods1.getQuick(j, i));
				// and the lower g rows with the a-periods of the dual mesh
				M.setQuick(i + g, j, aPeriods2.getQuick(j, i));
			}
		}

		// print(M, 4);

		return M;
	}

	/**
	 * Calculates 2*g harmonic differentials on on a Delaunay triangulated
	 * surface hds corresponding to a specified homology basis. The weight
	 * Adapter is used to find a basis of the homology that are short with
	 * respect to the weight given by wa. One differential corresponds to the
	 * entries in a row.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param homologyBasis
	 * @param adapters
	 * @param la
	 * @param wa
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix2D getHarmonicFormsOfPrimalMesh(
		HalfEdgeDataStructure<V, E, F> delaunay,
		List<List<E>> homologyBasis, 
		AdapterSet adapters
	) {
		// For each cycle in the basis we get a corresponding harmonic
		// differential (closed but not exact, its integral differs by one on
		// the cycle).
		DoubleMatrix2D dh = DoubleFactory2D.sparse.make(homologyBasis.size(), delaunay.numEdges() / 2);
		DoubleMatrix1D form;
		for (int i = 0; i < homologyBasis.size(); i++) {
			form = getStandardHarmonicFormOnPrimalMesh(delaunay, adapters,
					homologyBasis.get(i));
			for (int j = 0; j < form.size(); j++) {
				if (form.getQuick(j) != 0)
					dh.setQuick(i, j, form.getQuick(j));
			}
		}

		return dh;
	}

	/**
	 * Calculates 2*g harmonic differentials on on the dual a Delaunay
	 * triangulated surface hds corresponding to a specified homology basis. The
	 * weight Adapter is used to find a basis of the homology that are short
	 * with respect to the weight given by wa. One differential corresponds to
	 * the entries in a row.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param dualHomologyBasis
	 * @param adapters
	 * @param la
	 * @param wa
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix2D getHarmonicFormsOfDualMesh(
		HalfEdgeDataStructure<V, E, F> delaunay,
		List<List<E>> dualHomologyBasis, 
		AdapterSet adapters
	) {
		// For each cycle in the basis we get a corresponding harmonic
		// differential (closed but not exact, its integral differs by one on
		// the cycle).
		DoubleMatrix2D dh = DoubleFactory2D.sparse.make(dualHomologyBasis
				.size(), delaunay.numEdges() / 2);

		DoubleMatrix1D form;
		for (int i = 0; i < dualHomologyBasis.size(); i++) {
			form = getStandardHarmonicFormOnDualMesh(delaunay, adapters,
					dualHomologyBasis.get(i));
			for (int j = 0; j < form.size(); j++) {
				if (form.getQuick(j) != 0)
					dh.setQuick(i, j, form.getQuick(j));
			}
		}
		return dh;
	}

	/**
	 * Returns the standard harmonic differential on the primal mesh for the
	 * given cycle, i.e. the harmonic differential, whose integral along the
	 * cycle is 1.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param cycle
	 * @return
	 */
	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> DoubleMatrix1D getStandardHarmonicFormOnPrimalMesh(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle) {

		// get dual cycle, this is the chain of all edges ending at the left of
		// the given cycle
		List<E> edgesEndingAtLeft = DualityUtility.getDualPath(hds, cycle);

		// Consider the surface hds being obtained by identifying two boundary
		// cycles of another surface M via the map above. On this surface we can
		// get a harmonic function with boundary conditions 0 and 1 constant on
		// the identified cycles.
		DoubleMatrix1D h = getStandardHarmonicFunctionOnPrimalMesh(hds,
				adapters, cycle, edgesEndingAtLeft);

		// The differential of h can be defined on hds.
		DoubleMatrix1D dh = DoubleFactory1D.dense.make(hds.numEdges() / 2);

		// build the differences for each edge
		for (E e : hds.getPositiveEdges()) {
			int k = adapters.get(EdgeIndex.class, e, Integer.class);
			if (edgesEndingAtLeft.contains(e.getOppositeEdge())) {
				dh.setQuick(k, h.getQuick(e.getTargetVertex().getIndex()) + 1
						- h.getQuick(e.getStartVertex().getIndex()));
			} else if (edgesEndingAtLeft.contains(e)) {
				dh.setQuick(k, h.getQuick(e.getTargetVertex().getIndex())
						- h.getQuick(e.getStartVertex().getIndex()) - 1);
			} else {
				dh.setQuick(k, h.getQuick(e.getTargetVertex().getIndex())
						- h.getQuick(e.getStartVertex().getIndex()));
			}
		}

		System.out.println("Harmonic error: "
				+ howHarmonicIsPrimalForm(hds, adapters, dh));
		System.out.println();

		return dh;
	}

	/**
	 * Returns the standard harmonic differential on the dual mesh for the given
	 * cycle, i.e. the harmonic differential, whose integral along the cycle is
	 * 1.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param cycle
	 * @return
	 */
	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> DoubleMatrix1D getStandardHarmonicFormOnDualMesh(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle) {

		// get dual cycle, this is the chain of all edges ending at the left of
		// the given cycle
		List<E> edgesEndingAtLeft = DualityUtility.getPrimalPath(hds, cycle);

		// Consider the surface hds being obtained by identifying two boundary
		// cycles of another surface M via the map above. On this surface we can
		// get a harmonic function with boundary conditions 0 and 1 constant on
		// the identified cycles.
		DoubleMatrix1D h = getStandardHarmonicFunctionOnDualMesh(hds, adapters,
				cycle, edgesEndingAtLeft);

		// The differential of h can be defined on hds.
		DoubleMatrix1D dh = DoubleFactory1D.dense.make(hds.numEdges() / 2);

		// build the differences for each edge
		for (E e : hds.getPositiveEdges()) {
			int k = adapters.get(EdgeIndex.class, e, Integer.class);
			if (edgesEndingAtLeft.contains(e.getOppositeEdge())) {
				dh.setQuick(k, h.getQuick(e.getLeftFace().getIndex()) + 1
						- h.getQuick(e.getRightFace().getIndex()));
			} else if (edgesEndingAtLeft.contains(e)) {
				dh.setQuick(k, h.getQuick(e.getLeftFace().getIndex())
						- h.getQuick(e.getRightFace().getIndex()) - 1);
			} else {
				dh.setQuick(k, h.getQuick(e.getLeftFace().getIndex())
						- h.getQuick(e.getRightFace().getIndex()));
			}
		}

		System.out.println("Harmonic error: "
				+ howHarmonicIsDualForm(hds, adapters, dh));
		System.out.println();

		return dh;
	}

	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> double howHarmonicIsPrimalForm(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			DoubleMatrix1D form) {

		int id;
		double weight, res = 0;
		CotanWeightAdapter ca = new CotanWeightAdapter();
		adapters.add(ca);
		// int counter=0;
		for (V v : hds.getVertices()) {
			List<E> star = HalfEdgeUtilsExtra.getEdgeStar(v);
			double curr = 0;
			for (E e : star) {
				weight = adapters.get(Weight.class, e, Double.class);
				id = adapters.get(EdgeIndex.class, e, Integer.class);
				if (e.isPositive())
					curr += weight * form.getQuick(id);
				else
					curr -= weight * form.getQuick(id);
			}
			// if (Math.abs(curr) > 0.001)
			// System.out.println((counter++)+" Vertex " + v.getIndex() +
			// ": Value = "
			// + curr);
			res += Math.abs(curr);
		}
		adapters.remove(ca);
		return res;
	}

	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> double howHarmonicIsDualForm(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			DoubleMatrix1D form) {

		int id;
		double weight, res = 0;
		CotanWeightAdapter ca = new CotanWeightAdapter();
		adapters.add(ca);
		// int counter=0;
		for (F f : hds.getFaces()) {
			List<E> star = HalfEdgeUtilsExtra.getBoundary(f);
			double curr = 0;
			for (E e : star) {
				weight = 1. / adapters.get(Weight.class, e, Double.class);
				id = adapters.get(EdgeIndex.class, e, Integer.class);
				if (e.isPositive())
					curr -= weight * form.getQuick(id);
				else
					curr += weight * form.getQuick(id);
			}
			// if (Math.abs(curr) > 0.001)
			// System.out.println((counter++)+" Face " + f.getIndex() +
			// ": Value = "
			// + curr);
			res += Math.abs(curr);
		}
		adapters.remove(ca);
		return res;
	}

	/**
	 * Returns a harmonic function on the primal mesh which has a jump of 1 on
	 * the cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param cycle
	 * @param edgesEndingAtLeftOfCycle
	 * @return
	 */
	private static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> DoubleMatrix1D getStandardHarmonicFunctionOnPrimalMesh(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle, List<E> edgesEndingAtLeftOfCycle) {

		// the dimension of the linear system
		int n = hds.numVertices();

		DoubleMatrix2D laplaceop = LaplaceUtility.getPrimalLaplacian(hds,
				adapters);

		DoubleMatrix1D diag = DoubleFactory1D.dense.make(n);
		for (int i = 0; i < n; i++) {
			diag.setQuick(i, laplaceop.getQuick(i, i));
		}

		DoubleMatrix1D x = DoubleFactory1D.dense.make(n);

		double[] bcond = new double[n];
		double weight;

		CotanWeightAdapter ca = new CotanWeightAdapter();
		adapters.add(ca);

		// the function shall have a jump of 1 crossing the cycle
		for (E e : edgesEndingAtLeftOfCycle) {
			weight = adapters.get(Weight.class, e, Double.class);
			bcond[e.getStartVertex().getIndex()] += weight;
			bcond[e.getTargetVertex().getIndex()] -= weight;
		}

		adapters.remove(ca);

		DoubleMatrix1D b = DoubleFactory1D.dense.make(bcond);

		int J = -1;
		for (int i = 0; i < n; i++) {
			if (diag.getQuick(i) != 0) {
				for (int j = 0; j < n; j++) {
					if (laplaceop.getQuick(i, j) != 0)
						laplaceop.setQuick(i, j, laplaceop.getQuick(i, j)
								/ diag.getQuick(i));
				}
				b.setQuick(i, b.getQuick(i) / diag.getQuick(i));
			}
			if (J < 0 && b.getQuick(i) == 0)
				J = i;
		}

		// since the function is only unique up to constants we can fix the 0th
		// value to 1
		for (int i = 0; i < n; i++) {
			if (i == J)
				laplaceop.setQuick(J, i, 1);
			else
				laplaceop.setQuick(J, i, 0);
		}
		b.setQuick(J, 1);
		// System.err.println("determinant: "+ dalgebra.det(laplaceop));

		solve(laplaceop, x, b);

		DoubleMatrix1D H = DoubleFactory1D.dense.make(n);

		for (int i = 0; i < n; i++) {
			H.setQuick(i, x.getQuick(i));
		}

		return H;

	}

	/**
	 * Returns a harmonic function on the dual mesh which has a jump of 1 on the
	 * cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param cycle
	 * @param edgesEndingAtLeftOfCycle
	 * @return
	 */
	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> DoubleMatrix1D getStandardHarmonicFunctionOnDualMesh(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle, List<E> edgesEndingAtLeftOfCycle) {

		int n = hds.numFaces();

		DoubleMatrix2D laplaceop = LaplaceUtility.getDualLaplacian(hds,
				adapters);

		DoubleMatrix1D diag = DoubleFactory1D.dense.make(n);
		for (int i = 0; i < n; i++) {
			diag.setQuick(i, laplaceop.getQuick(i, i));
		}

		DoubleMatrix1D x = DoubleFactory1D.dense.make(n);

		double[] bcond = new double[n];
		double weight;

		CotanWeightAdapter ca = new CotanWeightAdapter();
		adapters.add(ca);

		// the function shall have a jump of 1 crossing the cycle
		for (E e : edgesEndingAtLeftOfCycle) {
			weight = 1. / adapters.get(Weight.class, e, Double.class);
			bcond[e.getRightFace().getIndex()] += weight;
			bcond[e.getLeftFace().getIndex()] -= weight;
		}

		adapters.remove(ca);

		DoubleMatrix1D b = DoubleFactory1D.dense.make(bcond);

		int J = -1;
		for (int i = 0; i < n; i++) {
			if (diag.getQuick(i) != 0) {
				for (int j = 0; j < n; j++) {
					if (laplaceop.getQuick(i, j) != 0)
						laplaceop.setQuick(i, j, laplaceop.getQuick(i, j)
								/ diag.getQuick(i));
				}
				b.setQuick(i, b.getQuick(i) / diag.getQuick(i));
				if (J < 0 && b.getQuick(i) == 0)
					J = i;
			}
		}

		// since the function is only unique up to constants we can fix the 0th
		// value to 1
		for (int i = 0; i < n; i++) {
			if (i == J)
				laplaceop.setQuick(J, i, 1);
			else
				laplaceop.setQuick(J, i, 0);
		}
		b.setQuick(J, 1);

		// System.err.println("determinant: "+ dalgebra.det(laplaceop));

		solve(laplaceop, x, b);

		DoubleMatrix1D H = DoubleFactory1D.dense.make(n);

		for (int i = 0; i < n; i++) {
			H.setQuick(i, x.getQuick(i));
		}

		return H;
	}

	// private static double eps = 1E-20;
	// private static int maxIterations = 100000000;

	/**
	 * Solves Ax=b and writes the result in the vector x.
	 * 
	 * @param A
	 * @param x
	 * @param b
	 */
	private static void solve(DoubleMatrix2D A, DoubleMatrix1D x,
			DoubleMatrix1D b) {

		// DoubleIterativeSolver solver;
		// // solver = new DoubleGMRES(x);
		// solver = new DoubleBiCGstab(x);
		//
		// DefaultDoubleIterationMonitor monitor = new
		// DefaultDoubleIterationMonitor();
		//
		// // configure monitor
		// monitor.setMaxIterations(maxIterations);
		// monitor.setAbsoluteTolerance(eps);
		// monitor.setRelativeTolerance(eps);
		// // monitor.setDivergenceTolerance(1);
		// monitor.setNormType(Norm.Two);
		// monitor.setIterationReporter(reporter);
		//
		// solver.setIterationMonitor(monitor);
		//
		// try {
		// solver.solve(A, b, x);
		// } catch (IterativeSolverDoubleNotConvergedException e) {
		// System.err
		// .println("Iterative solver failed to converge: Couldn't get harmonic function.");
		// e.printStackTrace();
		// }
		// System.err.println();

		System.out.println("Start solving");
		long startingtime = System.currentTimeMillis();
		
		for (int i = 0; i < b.size(); i++) {
			x.setQuick(i, b.getQuick(i));
		}

		SparseDoubleLUDecomposition lu = new SparseDoubleLUDecomposition(
				((SparseDoubleMatrix2D) A).getColumnCompressed(false), 1, false);
		lu.solve(x);
		
		System.out.println("time: "+ (System.currentTimeMillis()-startingtime));
	}

}
