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
import cern.colt.matrix.tdouble.algo.solver.DoubleIterationReporter;
import cern.colt.matrix.tdouble.algo.solver.DoubleIterativeSolver;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.EdgeIndex;
import de.jtem.halfedgetools.adapter.type.Weight;
import de.jtem.halfedgetools.algorithm.triangulation.Delaunay;
import de.jtem.halfedgetools.algorithm.triangulation.MappedLengthAdapter;
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

	private static DoubleIterationReporter reporter = new ColtIterationReporterImpl();	
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
		// First make clear that we are working with a delaunay triangulation.
		MappedLengthAdapter la = Delaunay.constructDelaunay(hds, adapters);

		// Get the homology basis of the surface.
		V rootV = hds.getVertex(0);
		List<List<E>> basis = CanonicalBasisUtility.getCanonicalHomologyBasis(
				rootV, adapters, wa);

		DoubleMatrix2D dh= getHarmonicFormsOfPrimalMesh(hds, basis, adapters, la);

		SimpleMatrixPrintUtility.print(dalgebra.mult(dh, EdgeUtility
				.cyclesToMatrix(adapters, hds, basis)), 4);

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
		// First make clear that we are working with a delaunay triangulation.
		MappedLengthAdapter la = Delaunay.constructDelaunay(hds, adapters);

		// Get the homology basis of the surface.
		V rootV = hds.getVertex(0);
		List<List<E>> dualbasis = DualityUtility.getDualPaths(hds,
				CanonicalBasisUtility.getCanonicalHomologyBasis(rootV,
						adapters, wa));

		DoubleMatrix2D dh = getHarmonicFormsOfDualMesh(hds, dualbasis,
				adapters, la);

		SimpleMatrixPrintUtility.print(dalgebra.mult(dh, EdgeUtility
				.cyclesToMatrix(adapters, hds, dualbasis)), 4);

		// use the private method
		return dh.toArray();
	}

	/**
	 * Returns a matrix A needed to build the normalized holomorphic
	 * differentials. The matrix has format 2g*2g. There are 2g harmonic
	 * differentials dH_j= dh1_j+i*dh2_j. The first g rows are filled with
	 * its a1-periods, the second g rows with its values on the a2-cycles. 
	 * That means in the upper g rows of Ax stands the real part and in
	 * the lower the imaginary part of the a-periods of the holomorphic
	 * differential omega=sum(x_j*dH_j).
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
	public static DoubleMatrix2D getHarmonicPeriodsMatrix(
			DoubleMatrix2D DH1, DoubleMatrix2D DH2, DoubleMatrix2D A1,
			DoubleMatrix2D A2) {

		// The genus equals the number of a-cycles equals the number of columns of a
		int g = A1.columns();

		// get a-period-matrix of the harmonic differentials DH
		DoubleMatrix2D aPeriods1 = dalgebra.mult(DH1,A1);
		DoubleMatrix2D aPeriods2 = dalgebra.mult(DH2,A2);

		// The matrix has to have the format 2g*2g. 
		DoubleMatrix2D M = DoubleFactory2D.sparse.make(2 * g, 2 * g);

		// for each cycle and each differential
		for (int i = 0; i < g; i++) {
			for (int j = 0; j < 2 * g; j++) {
				// fill the upper g rows with a-periods of the primal mesh
				M.set(i, j, aPeriods1.get(j, i));
				// and the lower g rows with the a-periods of the dual mesh
				M.set(i + g, j, aPeriods2.get(j, i));
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
				List<List<E>> homologyBasis, AdapterSet adapters,
				MappedLengthAdapter la) {
		
		// For each cycle in the basis we get a corresponding harmonic
		// differential (closed but not exact, its integral differs by one on
		// the cycle).
		DoubleMatrix2D dh = DoubleFactory2D.sparse.make(homologyBasis.size(),
				delaunay.numEdges() / 2);

		// For each cycle construct the corresponding harmonic differential.
		adapters.add(la); // delaunay lengths
		
		DoubleMatrix1D form;
		for (int i = 0; i < homologyBasis.size(); i++) {
			form= getStandardHarmonicFormOnPrimalMesh(delaunay, adapters,
					homologyBasis.get(i));
			for (int j = 0; j < form.size(); j++) {
				if (form.get(j) != 0)
					dh.set(i, j, form.get(j));
			}
		}

		adapters.remove(la);
		
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
				List<List<E>> dualHomologyBasis, AdapterSet adapters,
				MappedLengthAdapter la) {
		
		// For each cycle in the basis we get a corresponding harmonic
		// differential (closed but not exact, its integral differs by one on
		// the cycle).
		DoubleMatrix2D dh = DoubleFactory2D.sparse.make(dualHomologyBasis.size(),
				delaunay.numEdges() / 2);

		// For each cycle construct the corresponding harmonic differential.
		adapters.add(la); // delaunay lengths
		
		DoubleMatrix1D form;
		for (int i = 0; i < dualHomologyBasis.size(); i++) {
			form= getStandardHarmonicFormOnDualMesh(delaunay, adapters,
					dualHomologyBasis.get(i));
			for (int j = 0; j < form.size(); j++) {
				if (form.get(j) != 0)
					dh.set(i, j, form.get(j));
			}
		}

		adapters.remove(la);

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
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix1D getStandardHarmonicFormOnPrimalMesh(
				HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
				List<E> cycle) {

		// get dual cycle, this is the chain of all edges ending at the left of the given cycle
		List<E> edgesEndingAtLeft= DualityUtility.getDualPath(hds, cycle);
		
		// Consider the surface hds being obtained by identifying two boundary
		// cycles of another surface M via the map above. On this surface we can
		// get a harmonic function with boundary conditions 0 and 1 constant on
		// the identified cycles.
		DoubleMatrix1D h = getStandardHarmonicFunctionOnPrimalMesh(hds, adapters, cycle, edgesEndingAtLeft);

		// The differential of h can be defined on hds.
		DoubleMatrix1D dh = DoubleFactory1D.dense.make(hds.numEdges() / 2);

		// build the differences for each edge
		for (E e : hds.getPositiveEdges()) {
			int k = adapters.get(EdgeIndex.class, e, Integer.class);
			if (edgesEndingAtLeft.contains(e.getOppositeEdge())) {
				dh.set(k, h.get(e.getTargetVertex().getIndex()) + 1
						- h.get(e.getStartVertex().getIndex()));
			} else if (edgesEndingAtLeft.contains(e)) {
				dh.set(k, h.get(e.getTargetVertex().getIndex())
						- h.get(e.getStartVertex().getIndex()) - 1);
			} else {
				dh.set(k, h.get(e.getTargetVertex().getIndex())
						- h.get(e.getStartVertex().getIndex()));
			}
		}
		
		System.out.println("Harmonic error: "+howHarmonicIsPrimalForm(hds, adapters, dh));
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
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix1D getStandardHarmonicFormOnDualMesh(
				HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
				List<E> cycle) {

		// get dual cycle, this is the chain of all edges ending at the left of the given cycle
		List<E> edgesEndingAtLeft= DualityUtility.getPrimalPath(hds, cycle);
		
		// Consider the surface hds being obtained by identifying two boundary
		// cycles of another surface M via the map above. On this surface we can
		// get a harmonic function with boundary conditions 0 and 1 constant on
		// the identified cycles.
		DoubleMatrix1D h = getStandardHarmonicFunctionOnDualMesh(hds, adapters,
				cycle,edgesEndingAtLeft);

		// The differential of h can be defined on hds.
		DoubleMatrix1D dh = DoubleFactory1D.dense.make(hds.numEdges() / 2);

		// build the differences for each edge
		for (E e : hds.getPositiveEdges()) {
			int k = adapters.get(EdgeIndex.class, e, Integer.class);
			if (edgesEndingAtLeft.contains(e.getOppositeEdge())) {
				dh.set(k, h.get(e.getLeftFace().getIndex()) + 1
						- h.get(e.getRightFace().getIndex()));
			} else if (edgesEndingAtLeft.contains(e)) {
				dh.set(k, h.get(e.getLeftFace().getIndex())
						- h.get(e.getRightFace().getIndex()) - 1);
			} else {
				dh.set(k, h.get(e.getLeftFace().getIndex())
						- h.get(e.getRightFace().getIndex()));
			}
		}
		
		System.out.println("Harmonic error: "+howHarmonicIsDualForm(hds, adapters, dh));
		System.out.println();

		return dh;
	}
	
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double howHarmonicIsPrimalForm(
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
				if(e.isPositive())
					curr += weight * form.get(id);
				else
					curr -= weight * form.get(id);
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

	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double howHarmonicIsDualForm(
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
				weight = 1./adapters.get(Weight.class, e, Double.class);
				id = adapters.get(EdgeIndex.class, e, Integer.class);
				if(e.isPositive())
					curr -= weight * form.get(id);
				else
					curr += weight * form.get(id);
			}
			// if (Math.abs(curr) > 0.001)
			// System.out.println((counter++)+" Face " + f.getIndex() +
			// ": Value = "
			//		+ curr);
			res += Math.abs(curr);
		}
		adapters.remove(ca);
		return res;
	}
	
	/**
	 * Returns a harmonic function on the primal mesh which has a jump of 1 on the cycle.
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
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix1D getStandardHarmonicFunctionOnPrimalMesh(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle, List<E> edgesEndingAtLeftOfCycle) {

		// the dimension of the linear system
		int n = hds.numVertices();

		DoubleMatrix2D laplaceop = LaplaceUtility.getPrimalLaplacian(hds,
				adapters);
		
		DoubleMatrix1D diag= DoubleFactory1D.dense.make(n);
		for (int i = 0; i < n; i++) {
			diag.set(i, laplaceop.get(i, i));
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
			if(diag.get(i)!=0){
				for (int j = 0; j < n; j++) {
					if(laplaceop.get(i, j)!=0)
						laplaceop.set(i, j, laplaceop.get(i, j)/diag.get(i));
				}
				b.set(i, b.get(i)/diag.get(i));
			}
			if(J<0&&b.get(i)==0)
				J=i;
		}
		
		// since the function is only unique up to constants we can fix the 0th value to 1
		for (int i = 0; i < n; i++) {
			if(i==J)
				laplaceop.set(J, i, 1);
			else
				laplaceop.set(J, i, 0);
		}
		b.set(J, 1);
		// System.err.println("determinant: "+ dalgebra.det(laplaceop));

		solve(laplaceop, x, b);

		DoubleMatrix1D H = DoubleFactory1D.dense.make(n);

		for (int i = 0; i < n; i++) {
			H.set(i, x.get(i));
		}

		return H;
		
	}

	/**
	 * Returns a harmonic function on the dual mesh which has a jump of 1 on the cycle.
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
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix1D getStandardHarmonicFunctionOnDualMesh(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle, List<E> edgesEndingAtLeftOfCycle) {
		
		int n = hds.numFaces();
		
		DoubleMatrix2D laplaceop = LaplaceUtility.getDualLaplacian(hds, adapters);
		
		DoubleMatrix1D diag= DoubleFactory1D.dense.make(n);
		for (int i = 0; i < n; i++) {
			diag.set(i, laplaceop.get(i, i));
		}

		DoubleMatrix1D x = DoubleFactory1D.dense.make(n);

		double[] bcond = new double[n];
		double weight;

		CotanWeightAdapter ca = new CotanWeightAdapter();
		adapters.add(ca);

		// the function shall have a jump of 1 crossing the cycle
		for (E e : edgesEndingAtLeftOfCycle) {
			weight = 1./adapters.get(Weight.class, e, Double.class);
			bcond[e.getRightFace().getIndex()] += weight;
			bcond[e.getLeftFace().getIndex()] -= weight;
		}

		adapters.remove(ca);

		DoubleMatrix1D b = DoubleFactory1D.dense.make(bcond);
		
		int J= -1;
		for (int i = 0; i < n; i++) {
			if(diag.get(i)!=0){
				for (int j = 0; j < n; j++) {
					if(laplaceop.get(i, j)!=0)
						laplaceop.set(i, j, laplaceop.get(i, j)/diag.get(i));
				}
				b.set(i, b.get(i)/diag.get(i));
				if(J<0&&b.get(i)==0)
					J=i;
			}
		}
		
		// since the function is only unique up to constants we can fix the 0th value to 1
		for (int i = 0; i < n; i++) {
			if(i==J)
				laplaceop.set(J, i, 1);
			else
				laplaceop.set(J, i, 0);
		}
		b.set(J, 1);
		
		// System.err.println("determinant: "+ dalgebra.det(laplaceop));
		
		solve(laplaceop, x, b);

		DoubleMatrix1D H = DoubleFactory1D.dense.make(n);

		for (int i = 0; i < n; i++) {
			H.set(i, x.get(i));
		}

		return H;
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
		monitor.setNormType(Norm.Two);
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
