package de.varylab.discreteconformal.util;

import java.text.NumberFormat;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.Vector;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.Norm;
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.solver.DefaultDoubleIterationMonitor;
import cern.colt.matrix.tdouble.algo.solver.DoubleBiCG;
import cern.colt.matrix.tdouble.algo.solver.DoubleBiCGstab;
import cern.colt.matrix.tdouble.algo.solver.DoubleGMRES;
import cern.colt.matrix.tdouble.algo.solver.DoubleIterationReporter;
import cern.colt.matrix.tdouble.algo.solver.DoubleIterativeSolver;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.EdgeIndex;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.Weight;
import de.jtem.halfedgetools.algorithm.triangulation.Delaunay;
import de.jtem.halfedgetools.algorithm.triangulation.MappedLengthAdapter;
import de.jtem.halfedgetools.util.HalfEdgeUtilsExtra;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.util.Search.WeightAdapter;


/**
 * Class to calculate harmonic differentials corresponding to a cycle of the
 * homology basis and holomorphic differentials in the sense of mercat.
 * 
 * By convention all the forms are saved as rows and all the cycles as columns.
 * (For private methods we always use colt matrices.)
 * 
 * @author knoeppel
 * 
 */
public class DiscreteRiemannUtility {

	/**
	 * Adapter returns the cotan weight of a given edge.
	 */
	@Weight
	private static class CotanAdapter extends AbstractAdapter<Double> {

		public CotanAdapter() {
			super(Double.class, true, false);
		}

		public <
			V extends Vertex<V, E, F>,
			E extends Edge<V, E, F>,
			F extends Face<V, E, F>
		> Double getE(E e, AdapterSet adapters) {
			if (!adapters.contains(Length.class, e.getClass(), Double.class)) {
				throw new RuntimeException(
						"Need adapter for length of edges to calculate cotan weights.");
			}

			double a = adapters.get(Length.class, e, Double.class);
			double b = adapters
					.get(Length.class, e.getNextEdge(), Double.class);
			double c = adapters.get(Length.class,
					e.getNextEdge().getNextEdge(), Double.class);

			if (!e.getNextEdge().getNextEdge().getNextEdge().equals(e)) {
				throw new RuntimeException("Face is not a triangle.");
			}

			double cosalpha = (b * b + c * c - a * a) / (2 * b * c);
			double cotanalpha = cosalpha / (Math.sqrt(1 - cosalpha * cosalpha));

			return cotanalpha;
		};

		@Override
		public <
			N extends Node<?, ?, ?>
		> boolean canAccept(Class<N> nodeClass) {
			return Edge.class.isAssignableFrom(nodeClass);
		}

	}

	private static DoubleIterationReporter reporter = new DoubleIterationReporter() {
		
		@Override
		public void monitor(double arg0, DoubleMatrix1D arg1, int arg2) {
			monitor(arg0, arg2);
		}

		@Override
		public void monitor(double arg0, int arg1) {
			if (arg1 % 100 == 0)
				System.err.println("iteration = " + arg1 + ", value = " + arg0);
		}
	};
	
	private static DenseDoubleAlgebra dalgebra = new DenseDoubleAlgebra();
	private static double eps = 1E-10;
	
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
		List<List<E>> basis = getCanonicalHomologyBasis(rootV,adapters, wa);
		
		DoubleMatrix2D dh= getHarmonicFormsOfPrimalMesh(hds, basis, adapters, la, wa);

		// print(dalgebra.mult(dh, cyclesToMatrix(adapters, hds,basis)), 4);
		
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
		List<List<E>> dualbasis = getDualPaths(hds, getCanonicalHomologyBasis(rootV,adapters, wa));
		
		DoubleMatrix2D dh= getHarmonicFormsOfDualMesh(hds, dualbasis, adapters, la, wa);

		// print(dalgebra.mult(dh, cyclesToMatrix(adapters, hds,basis)), 4);
		
		// use the private method
		return dh.toArray();
	}

	/**
	 * Returns a basis of holomorphic differentials on the surface hds. The rows
	 * represent holomorphic forms. In the i-th column the value on the positive
	 * oriented edge with edgeIndex i is saved. The real part correspond to the
	 * edge the imaginary to its dual.
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
	> Complex[][] getHolomorphicForms(
		HalfEdgeDataStructure<V, E, F> hds, 
		AdapterSet adapters,
		WeightAdapter<E> wa){
		
		// First make clear that we are working with a delaunay triangulation.
		MappedLengthAdapter la = Delaunay.constructDelaunay(hds, adapters);

		// Get the homology basis of the surface.
		V rootV = hds.getVertex(0);
		List<List<E>> basis = getCanonicalHomologyBasis(rootV,adapters, wa);
		List<List<E>> dualbasis = getDualPaths(hds,basis);
		
		DoubleMatrix2D[] omega1= getHolomorphicFormsOnPrimalMesh(hds, basis, adapters, la, wa);
		DoubleMatrix2D[] omega2= getHolomorphicFormsOnDualMesh(hds, dualbasis, adapters, la, wa);
		
		// to normalize the differentials we need the a-periods and its dual
		// cycles
		List<List<E>> acycles = getACycles(basis);
		List<List<E>> dualacycles = getACycles(dualbasis);

		// write cycles to matrices
		DoubleMatrix2D A = cyclesToMatrix(adapters, hds,acycles);
		DoubleMatrix2D dualA = cyclesToMatrix(adapters, hds, dualacycles);

		System.err.println();
		System.err.println("A - PERIODS:");
		System.out.println("real part of primal construction:");
		print(dalgebra.mult(omega1[0], A), 4);
		System.err.println();
		System.out.println("imaginary part of primal construction:");
		print(dalgebra.mult(omega1[1], dualA), 4);
		System.err.println();
		System.err.println("A - PERIODS:");
		System.out.println("real part of dual construction:");
		print(dalgebra.mult(omega2[0], dualA), 4);
		System.err.println();
		System.out.println("imaginary part of dual construction:");
		print(dalgebra.mult(omega2[1], A), 4);
		System.err.println();
		
		// to normalize the differentials we need the a-periods and its dual
		// cycles
		List<List<E>> bcycles = getBCycles(basis);
		List<List<E>> dualbcycles = getDualPaths(hds, bcycles);

		// write cycles to matrices
		DoubleMatrix2D B = cyclesToMatrix(adapters, hds,bcycles);
		DoubleMatrix2D dualB = cyclesToMatrix(adapters, hds, dualbcycles);

		System.err.println();
		System.err.println("PERIOD MATRIX:");
		System.out.println("real part of primal construction:");
		print(dalgebra.mult(omega1[0], B), 4);
		System.err.println();
		System.out.println("imaginary part of primal construction:");
		print(dalgebra.mult(omega1[1], dualB), 4);
		System.err.println();
		System.err.println("PERIOD MATRIX:");
		System.out.println("real part of dual construction:");
		print(dalgebra.mult(omega2[0], dualB), 4);
		System.err.println();
		System.out.println("imaginary part of dual construction:");
		print(dalgebra.mult(omega2[1], B), 4);
		System.err.println();
		
		int m= omega1[0].rows(); int n= omega1[0].columns();
		Complex[][] array= new Complex[m][n];
		// TODO: check the minus sign this after all methods are doubled
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				array[i][j] = new Complex(omega1[0].get(i, j)
						- omega2[1].get(i, j), omega1[1].get(i, j)
						+ omega2[0].get(i, j)); 
			}
		}
		
		// use the private method
		return array; 
	}
	
	/**
	 * class representing a bilinear form
	 */
	private static class BilinearForm {
		private DoubleMatrix2D S;
		private int dim;

		public BilinearForm(DoubleMatrix2D S) {
			setMatrix(S);
		}
		
		public int getDimension() {
			return dim;
		}

		public void setMatrix(DoubleMatrix2D S) {
			// matrix must be square
			if (S.rows() != S.columns())
				throw new IllegalArgumentException(
						"Bilinear form couldn't be build: the matrix must be square.");
			this.dim = S.rows();
			this.S = S;
		}

		public double brackets(DoubleMatrix1D x, DoubleMatrix1D y) {
			// vector lengths must fit
			if (x.size() != dim || y.size() != dim)
				throw new IllegalArgumentException(
						"Couldn't evaluate brackets: vector dimension doesn't fit.");
			return dalgebra.mult(x, dalgebra.mult(S, y));
		}
	}

	/**
	 * Calculates a canonical homology basis {a_1,...,a_g,b_1,...,b_g}.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param root
	 * @param wa
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<List<E>> getCanonicalHomologyBasis(V root, AdapterSet adapters, WeightAdapter<E> wa){

		// get a homology basis (not canonical)
		List<Set<E>> paths = HomologyUtility.getGeneratorPaths(root, wa);

		// write the sets as cycles, i.e. lists (an edge can appear more than
		// only once)
		List<List<E>> homologyBasis = new Vector<List<E>>();
		for (Set<E> path : paths) {
			Vector<E> cycle = new Vector<E>();
			for (E e : path) {
				cycle.add(e);
			}
			homologyBasis.add(cycle);
		}

		// dimension of the homology group
		int dimension = homologyBasis.size();

		DoubleMatrix2D S = DoubleFactory2D.sparse.make(dimension, dimension);

		for (int i = 0; i < dimension; i++) {
			for (int j = i + 1; j < dimension; j++) {
				int s = getIntersectionNumberOfPrimalCycles(homologyBasis.get(i),
						homologyBasis.get(j));
				if (s != 0) {
					S.set(i, j, s);
					S.set(j, i, -s);
				}
			}
		}

		// print(S, 0);
		 
		DoubleMatrix2D Basis= cyclesToMatrix(adapters, root.getHalfEdgeDataStructure(), homologyBasis);
		
		// get coordinates of a canonical basis
		DoubleMatrix2D coords = canonicalBasisCoords(new BilinearForm(S));
		
		// build up new homology basis
		
		// calculate basis vectors
		DoubleMatrix2D BPrime = dalgebra.mult(Basis, coords);

		// convert to list representation
		List<List<E>> newHomologyBasis = matrixToCycles(adapters, root.getHalfEdgeDataStructure(), BPrime);

		// // TEST: show intersection number for canonical basis
		// System.out
		// .println("Intersection matrix for canonical homology basis:");
		// S = DoubleFactory2D.sparse.make(dimension, dimension);
		// for (int i = 0; i < dimension; i++) {
		// for (int j = i + 1; j < dimension; j++) {
		// int s = getIntersectionNumber(newHomologyBasis.get(i),
		// newHomologyBasis.get(j));
		// if (s != 0) {
		// S.set(i, j, s);
		// S.set(j, i, -s);
		// }
		// }
		// }
		//
		// print(S, 0);
		
		return newHomologyBasis;
		
	}
	
	/**
	 * Calculates the coordinates of a canonical basis of homology to a given
	 * basis with intersection form B.
	 * 
	 * @param B
	 * @return
	 */
	private static DoubleMatrix2D canonicalBasisCoords(BilinearForm B) {
		
		// dimension of the homology group
		int dimension = B.getDimension();

		// create a coordinate matrix
		DoubleMatrix2D coords = DoubleFactory2D.dense
				.make(dimension, dimension);
		for (int i = 0; i < dimension; i++) {
			coords.set(i, i, 1);
		}

		// initialize empty index sets
		Stack<Integer> todo = new Stack<Integer>();
		for (int i = 0; i < dimension; i++) {
			todo.add(i);
		}

		DoubleMatrix2D T;

		while (!todo.isEmpty()) {
			// get the smallest projection != 0
			Integer I = todo.pop();
			Integer J = getIndexOfSmallest(B, coords, I);

			if (J == -1)
				break;

			// remove the second index from index stack
			todo.remove(J);

			// build up basis transformation matrix
			T = DoubleFactory2D.sparse.make(dimension, dimension);
			for (int k = 0; k < dimension; k++) {
				if (k == I || k == J) {
					T.set(k, k, 1);
				} else {
					T.set(k,
							k,
							B.brackets(coords.viewColumn(I),
									coords.viewColumn(J)));
					T.set(J,
							k,
							B.brackets(coords.viewColumn(k),
									coords.viewColumn(I)));
					T.set(I,
							k,
							B.brackets(coords.viewColumn(J),
									coords.viewColumn(k)));
				}
			}
			coords = dalgebra.mult(coords, T);

			// System.err.println("new coords = ");
			// print(coords, 0);
			// System.err.println();
		}

		DoubleMatrix2D finalCoords = DoubleFactory2D.dense.make(dimension,
				dimension);

		DoubleMatrix2D S = DoubleFactory2D.sparse.make(dimension, dimension);
		for (int k = 0; k < dimension; k++) {
			for (int l = k + 1; l < dimension; l++) {
				double val = B.brackets(coords.viewColumn(k),
						coords.viewColumn(l));
				S.set(k, l, val);
				S.set(l, k, -val);
			}
		}
		
		// System.err.println("intersection matrix (not ordered)");
		// print(S, 0);

		IntArrayList rows = new IntArrayList();
		IntArrayList cols = new IntArrayList();
		DoubleArrayList vals = new DoubleArrayList();
		S.getPositiveValues(rows, cols, vals);

		for (int j = 0; j < rows.size(); j++) {
			for (int i = 0; i < dimension; i++) {
				finalCoords.set(i, j, coords.get(i, rows.get(j)));
				finalCoords.set(i, j + cols.size(), coords.get(i, cols.get(j)));
			}
		}

		// S= DoubleFactory2D.sparse.make(dimension,dimension);
		// for (int k = 0; k < dimension; k++) {
		// for (int l = k + 1; l < dimension; l++) {
		// double val = B.brackets(finalCoords.viewColumn(k),
		// finalCoords.viewColumn(l));
		// S.set(k, l, val);
		// S.set(l, k, -val);
		// }
		// }
		//
		// System.err.println("intersection matrix = ");
		// print(S,0);
	
		return finalCoords;
	}
	
	/**
	 * Returns the index of the cycle which has the smallest (absolute)
	 * intersection number != 0 with the cycle represented by the I-th column in the
	 * coordinate matrix.
	 */
	private static Integer getIndexOfSmallest(BilinearForm B,
			DoubleMatrix2D coords, Integer I) {
		Integer J = -1;
		double min = Double.MAX_VALUE;
		double curr;
		DoubleMatrix1D col = coords.viewColumn(I);
		for (int k = 0; k < B.getDimension(); k++) {
			curr = Math.abs(B.brackets(col, coords.viewColumn(k)));
			if (curr < min && curr != 0) {
				min = curr;
				J = k;
			}
		}
		return J;
	}

	/**
	 * Print method for colt matrices
	 * @param M
	 * @param noOfFractionDigits
	 */
	private static void print(DoubleMatrix2D M, int noOfFractionDigits) {
		
		NumberFormat numberFormat = NumberFormat
		.getInstance(Locale.ENGLISH);
		
		numberFormat.setMinimumFractionDigits(noOfFractionDigits);
		numberFormat.setMaximumFractionDigits(noOfFractionDigits);
		String tmp;
		for (int i = 0; i < M.rows(); i++) {
			System.out.print("| ");
			for (int j = 0; j < M.columns(); j++) {
				if (M.get(i, j) == 0)
					tmp = " " + numberFormat.format(0.);
				else {
					tmp = M.get(i, j) >= 0 ? "+" : "";
					tmp += numberFormat.format(M.get(i, j));
				}
				System.out.print(tmp);
				if (j < M.columns() - 1) {
					System.out.print(" , ");
				} else {
					System.out.println(" |");
				}
			}
		}
		System.out.println();
	}
	
	/**
	 * Print method for colt vector
	 * @param V
	 * @param noOfFractionDigits
	 */
	private static void print(DoubleMatrix1D V, int noOfFractionDigits) {

		NumberFormat numberFormat = NumberFormat.getInstance(Locale.ENGLISH);

		numberFormat.setMinimumFractionDigits(noOfFractionDigits);
		numberFormat.setMaximumFractionDigits(noOfFractionDigits);
		String tmp;

		System.out.print("| ");
		for (int i = 0; i < V.size(); i++) {
			if (V.get(i) == 0)
				tmp = " " + numberFormat.format(0.);
			else {
				tmp = V.get(i) >= 0 ? "+" : "";
				tmp += numberFormat.format(V.get(i));
			}
			System.out.print(tmp);
			if (i < V.size() - 1) {
				System.out.print(" , ");
			} else {
				System.out.println(" |");
			}
		}

		System.out.println();
	}

	/**
	 * Returns paths in the dual surface which are homotopic to the given
	 * ones.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param cycles
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<List<E>> getDualPaths(
				HalfEdgeDataStructure<V, E, F> hds, List<List<E>> cycles) {
		
		List<List<E>> dualCycles = new java.util.Vector<List<E>>();
		// for each cycle
		for (int i = 0; i < cycles.size(); i++) {
			dualCycles.add(getDualPath(hds, cycles.get(i)));
		}
		return dualCycles;
	}
	
	/**
	 * Returns a path in the dual surface which is homotopic to the given one.
	 * 	
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param cycle
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<E> getDualPath(
				HalfEdgeDataStructure<V, E, F> hds, List<E> cycle) {
		
		List<E> dualPath = new Vector<E>();
		// get vertices contained in the primal cycle
		Set<V> vertices = getPrimalVertexSet(cycle);
		// for each vertex in the cycle
		for (V v : vertices) {
			// get the edge star
			List<E> star = HalfEdgeUtilsExtra.getEdgeStar(v);
			// and test each edge in it
			for (E e : star) {
				EdgeStatus status = getPrimalEdgeStatus(e, cycle, vertices);
				// if the edge is on the left side put it to the dual path
				// and points to the vertex
				if (status == EdgeStatus.endsAtLeftCycle)
					dualPath.add(e);
				else if (status == EdgeStatus.startsAtLeftCycle)
					dualPath.add(e.getOppositeEdge());
				// dualPath.add(e);
			}
		}
		return dualPath;
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
	private static DoubleMatrix2D getHarmonicPeriodsMatrix(
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
	 * Gets canonical homology basis and returns g cycles of the homology basis
	 * representing the a-cycles, i.e. the first g cycles.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param canonicalHomologyBasis
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<List<E>> getACycles(List<List<E>> canonicalHomologyBasis){
		
		// the genus has to be one at least
		if (canonicalHomologyBasis.size() < 1)
			return null;

		// g non-intersecting cycles shall be labeled as a-cycles
		List<List<E>> aCycles = new Vector<List<E>>();
		
		for (int i= 0; i< canonicalHomologyBasis.size()/2; i++) {
			aCycles.add(i, canonicalHomologyBasis.get(i));
		}

		return aCycles;
	}
	
	/**
	 * Gets canonical homology basis and returns g cycles of the homology basis
	 * representing the b-cycles, i.e. the second g cycles.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param canonicalHomologyBasis
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<List<E>> getBCycles(List<List<E>> canonicalHomologyBasis){
		
		// the genus has to be one at least
		if (canonicalHomologyBasis.size() < 1)
			return null;

		// g non-intersecting cycles shall be labeled as a-cycles
		List<List<E>> bCycles = new Vector<List<E>>();
		
		for (int i= 0; i< canonicalHomologyBasis.size()/2; i++) {
			bCycles.add(i, canonicalHomologyBasis.get(i+canonicalHomologyBasis.size()/2));
		}

		return bCycles;
	}
	
	/**
	 * Returns the intersection number of two oriented primal cycles.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param cycle1
	 * @param cycle2
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> int getIntersectionNumberOfPrimalCycles(List<E> cycle1, List<E> cycle2) {
		
		// get the number of positive intersections of the first cycle with the
		// second
		int numPositiveIntersections = countPrimalEdgesWithStatus(cycle1, cycle2,
				EdgeStatus.startsAtLeftCycle);
		// get the number of negative intersections of the first cycle with the
		// second
		int numNegativeIntersections = countPrimalEdgesWithStatus(cycle1, cycle2,
				EdgeStatus.endsAtLeftCycle);
		return numPositiveIntersections - numNegativeIntersections;
	}
	
	/**
	 * Returns the intersection number of two oriented dual cycles.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param cycle1
	 * @param cycle2
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> int getIntersectionNumberOfDualCycles(List<E> cycle1, List<E> cycle2) {
		
		// get the number of positive intersections of the first cycle with the
		// second
		int numPositiveIntersections = countDualEdgesWithStatus(cycle1, cycle2,
				EdgeStatus.startsAtLeftCycle);
		// get the number of negative intersections of the first cycle with the
		// second
		int numNegativeIntersections = countDualEdgesWithStatus(cycle1, cycle2,
				EdgeStatus.endsAtLeftCycle);
		return numPositiveIntersections - numNegativeIntersections;
	}

	/**
	 * Count the edges in a given list having a specified status with respect to
	 * a given primal cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param cycle
	 * @param list
	 * @param status
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> int countPrimalEdgesWithStatus(
				List<E> cycle, List<E> list, EdgeStatus status) {
		
		int count = 0;
		Set<V> vertexSet= getPrimalVertexSet(cycle);
		for(E e:list){
			if(getPrimalEdgeStatus(e, cycle,vertexSet)== status)
				count++;
		}
		return count;
	}
	
	/**
	 * Count the edges in a given list having a specified status with respect to
	 * a given dual cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param cycle
	 * @param list
	 * @param status
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> int countDualEdgesWithStatus(
				List<E> cycle, List<E> list, EdgeStatus status) {
		
		int count = 0;
		Set<F> vertexSet= getDualVertexSet(cycle);
		for(E e:list){
			if(getDualEdgeStatus(e, cycle,vertexSet)== status)
				count++;
		}
		return count;
	}
	
	/**
	 * Applies the Hodge star operator to an array of forms on the surface and
	 * returns the result (forms are stored as rows).
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param adapters
	 * @param forms
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix2D getDualFromPrimalForms(
				HalfEdgeDataStructure<V, E, F> delaunay, AdapterSet adapters,
				DoubleMatrix2D forms) {
		
		DoubleMatrix2D formsStar = DoubleFactory2D.dense.make(forms.rows(), forms.columns());
		DoubleMatrix1D row;
		for (int i = 0; i < formsStar.rows(); i++) {
			row= getDualFromPrimalForm(delaunay, adapters, forms.viewRow(i));
			for (int j = 0; j < forms.columns(); j++) {
				formsStar.set(i, j, row.get(j));	
			}
		}
		
		return formsStar;
	}
	
	/**
	 * Returns how harmonic a differential is.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param form
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double howHarmonic(HalfEdgeDataStructure<V,E,F> hds, AdapterSet adapters, DoubleMatrix1D form){
		double abssum= 0;
		int id;
		for (V v: hds.getVertices()) {
			List<E> star= HalfEdgeUtilsExtra.getEdgeStar(v);
			double val= 0;
			for (E e: star) {
				id= adapters.get(EdgeIndex.class, e,Integer.class);
				if(e.isPositive())
					val+= form.get(id);
				else
					val-= form.get(id);
			}
			abssum+=Math.abs(val);
		}
		return abssum;
	}
	
	/**
	 * Applies the Hodge star operator to an array of forms on the dual surface and
	 * returns the result (forms are stored as rows).
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param adapters
	 * @param forms
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix2D getDualFromDualForms(
				HalfEdgeDataStructure<V, E, F> delaunay, AdapterSet adapters,
				DoubleMatrix2D forms) {
		
		DoubleMatrix2D formsStar = DoubleFactory2D.dense.make(forms.rows(), forms.columns());
		DoubleMatrix1D row;
		for (int i = 0; i < formsStar.rows(); i++) {
			row= getDualFromDualForm(delaunay, adapters, forms.viewRow(i));
			for (int j = 0; j < forms.columns(); j++) {
				formsStar.set(i, j, row.get(j));	
			}
		}
		
		return formsStar;
	}
	
	/**
	 * Applies the Hodge star operator to a form on the surface and
	 * returns the result.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param adapters
	 * @param form
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix1D getDualFromPrimalForm(
				HalfEdgeDataStructure<V, E, F> delaunay, AdapterSet adapters,
				DoubleMatrix1D form) {
		
		DoubleMatrix1D formStar = DoubleFactory1D.dense.make((int)form.size());

		// for the construction of the dual forms the cotan weights are needed
		CotanAdapter ca = new CotanAdapter();
		adapters.add(ca);
		
		double ratio;
		int id;

		for (E e : delaunay.getPositiveEdges()) {
			ratio = getCotanWeight(e, adapters);
			id = adapters.get(EdgeIndex.class, e, Integer.class);
			formStar.set(id, ratio * form.get(id));
		}

		adapters.remove(ca);
		
		return formStar;
	}
	
	/**
	 * Applies the Hodge star operator to a form on the dual surface and
	 * returns the result.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param adapters
	 * @param form
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix1D getDualFromDualForm(
				HalfEdgeDataStructure<V, E, F> delaunay, AdapterSet adapters,
				DoubleMatrix1D form) {
		
		DoubleMatrix1D formStar = DoubleFactory1D.dense.make((int)form.size());

		// for the construction of the dual forms the cotan weights are needed
		CotanAdapter ca = new CotanAdapter();
		adapters.add(ca);
		
		double ratio;
		int id;

		for (E e : delaunay.getPositiveEdges()) {
			ratio = getCotanWeight(e, adapters);
			id = adapters.get(EdgeIndex.class, e, Integer.class);
			// TODO: perhaps with minus?
			formStar.set(id, form.get(id) / ratio);
		}

		adapters.remove(ca);
		
		return formStar;
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
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix2D getHarmonicFormsOfPrimalMesh(
				HalfEdgeDataStructure<V, E, F> delaunay,
				List<List<E>> homologyBasis, AdapterSet adapters,
				MappedLengthAdapter la, WeightAdapter<E> wa) {
		
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
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix2D getHarmonicFormsOfDualMesh(
				HalfEdgeDataStructure<V, E, F> delaunay,
				List<List<E>> dualHomologyBasis, AdapterSet adapters,
				MappedLengthAdapter la, WeightAdapter<E> wa) {
		
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
		
		DoubleMatrix2D Cycles= cyclesToMatrix(adapters, delaunay, dualHomologyBasis);
		print(dalgebra.mult(dh, Cycles),4);
		
		adapters.remove(la);
		
		return dh;
	}

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
		DoubleMatrix2D dh = getHarmonicFormsOfPrimalMesh(delaunay, canonicalHomologyBasis, adapters, la,
				wa);
		DoubleMatrix2D dhStar = getDualFromPrimalForms(delaunay, adapters, dh);

		// to normalize the differentials we need the a-periods and its dual
		// cycles
		List<List<E>> acycles = getACycles(canonicalHomologyBasis);
		List<List<E>> dualacycles = getDualPaths(delaunay, acycles);

		// write cycles to matrices
		DoubleMatrix2D A = cyclesToMatrix(adapters, delaunay, acycles);
		DoubleMatrix2D dualA = cyclesToMatrix(adapters, delaunay, dualacycles);
		
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
		DoubleMatrix2D M = getHarmonicPeriodsMatrix(dh, dhStar,A,dualA);

		// print(M, 2);

		// The number of harmonic forms on the surface is 2g, so we need 2g
		// coefficients.
		DoubleMatrix1D x = DoubleFactory1D.dense.make(2 * g);

		// Iterative solver
		DoubleIterativeSolver solver = new DoubleGMRES(x);
		DefaultDoubleIterationMonitor monitor= new DefaultDoubleIterationMonitor();
		
		// configure monitor
		monitor.setMaxIterations(100000);
		monitor.setAbsoluteTolerance(eps);
		monitor.setRelativeTolerance(eps);
		monitor.setIterationReporter(reporter);
		
		monitor.setDivergenceTolerance(1);
		
		solver.setIterationMonitor(monitor);
				
		// For each a-cycle find the holomorphic form which is 1 along this
		// cycle, i.e. gives one for the primal cycle and 0 for its dual cycle. 
		for (int i = 0; i < g; i++) {
			
			// set up the conditions 
			DoubleMatrix1D bc = DoubleFactory1D.dense.make(2 * g);
			bc.set(g+i, 2*Math.PI);
			
			// solve the system
			try {
				solver.solve(M, bc, x);
			} catch (IterativeSolverDoubleNotConvergedException e) {
				System.err
						.println("Iterative solver failed to converge: Couldn't get holomorphic form.");
				e.printStackTrace();
			}

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

		// Get the harmonic differentials on the surface and its dual. The
		// format of the matrices is 2g*numEdges
		DoubleMatrix2D dhStar = getHarmonicFormsOfDualMesh(delaunay, canonicalHomologyBasis, adapters, la,
				wa);
		DoubleMatrix2D dhStarStar = getDualFromDualForms(delaunay, adapters, dhStar);

		// to normalize the differentials we need the a-periods and its dual
		// cycles
		List<List<E>> acycles = getACycles(canonicalHomologyBasis);
		List<List<E>> dualacycles = getDualPaths(delaunay, acycles);

		// write cycles to matrices
		DoubleMatrix2D A = cyclesToMatrix(adapters, delaunay, acycles);
		DoubleMatrix2D dualA = cyclesToMatrix(adapters, delaunay, dualacycles);
		
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
		DoubleMatrix2D M = getHarmonicPeriodsMatrix(dhStarStar, dhStar,dualA,A);

		// print(M, 2);

		// The number of harmonic forms on the surface is 2g, so we need 2g
		// coefficients.
		DoubleMatrix1D x = DoubleFactory1D.dense.make(2 * g);

		// Iterative solver
		DoubleIterativeSolver solver = new DoubleGMRES(x);
		DefaultDoubleIterationMonitor monitor= new DefaultDoubleIterationMonitor();
		
		// configure monitor
		monitor.setMaxIterations(100000);
		monitor.setAbsoluteTolerance(eps);
		monitor.setRelativeTolerance(eps);
		monitor.setIterationReporter(reporter);
		
		monitor.setDivergenceTolerance(1);
		
		solver.setIterationMonitor(monitor);
		
		// to normalize the differentials we need the a-periods and its dual
		// cycles
		List<List<E>> bcycles = getBCycles(canonicalHomologyBasis);
		List<List<E>> dualbcycles = getDualPaths(delaunay, acycles);

		// write cycles to matrices
		DoubleMatrix2D B = cyclesToMatrix(adapters, delaunay,bcycles);
		DoubleMatrix2D dualB = cyclesToMatrix(adapters, delaunay, dualbcycles);
		
		// For each a-cycle find the holomorphic form which is 1 along this
		// cycle, i.e. gives one for the primal cycle and 0 for its dual cycle. 
		for (int i = 0; i < g; i++) {
			
			// set up the conditions 
			DoubleMatrix1D bc = DoubleFactory1D.dense.make(2 * g);
			bc.set(g+i, 2*Math.PI);
			
			// solve the system
			try {
				solver.solve(M, bc, x);
			} catch (IterativeSolverDoubleNotConvergedException e) {
				System.err
						.println("Iterative solver failed to converge: Couldn't get holomorphic form.");
				e.printStackTrace();
			}

			// Build linear combinations of the harmonic differentials
			// using the obtained coefficients. 
			DoubleMatrix1D omega= dalgebra.mult(dalgebra.transpose(dhStar), x);
			DoubleMatrix1D omegaStar= dalgebra.mult(dalgebra.transpose(dhStarStar), x);
			
			for (int j = 0; j < numPosEdges; j++) {
				OMEGA[0].set(i, j, omega.get(j));
				OMEGA[1].set(i, j, omegaStar.get(j));
			}
		}
		
		System.err.println();
		System.err.println("PERIOD MATRIX:");
		System.out.println("real part:");
		print(dalgebra.mult(OMEGA[0], dualB), 4);
		System.err.println();
		System.out.println("imaginary part:");
		print(dalgebra.mult(OMEGA[1], B), 4);
		System.err.println();
		
		adapters.remove(la);

		return OMEGA;
	}

	/**
	 * Method gets a list of cycles and writes them into the columns of a matrix.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param adapters
	 * @param delaunay
	 * @param cycles
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix2D cyclesToMatrix(
			AdapterSet adapters, HalfEdgeDataStructure<V, E, F> delaunay,
			List<List<E>> cycles) {
		
		// initialize two matrices to encode the cycles
		DoubleMatrix2D A = new SparseDoubleMatrix2D(delaunay.numEdges() / 2,
				cycles.size());

		// needed to fill the matrices
		int currEdgeId;
		double currval;

		// fill the matrix of a-cycles
		for (int i = 0; i < cycles.size(); i++) {
			for (E e : cycles.get(i)) {
				currEdgeId = adapters.get(EdgeIndex.class, e, Integer.class);
				currval = A.get(currEdgeId, i);
				if (e.isPositive())
					currval++;
				else
					currval--;
				A.set(currEdgeId, i, currval);
			}
		}
		// print(A, 0);

		return A;
	}
	
	/**
	 * Method gets a matrix of cycles and returns them as a list.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param adapters
	 * @param delaunay
	 * @param cycles
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<List<E>> matrixToCycles(
			AdapterSet adapters, HalfEdgeDataStructure<V, E, F> delaunay,
			DoubleMatrix2D cycles) {
		
		// initialize two matrices to encode the cycles
		List<List<E>> A = new Vector<List<E>>();
		
		int n= cycles.columns();
		for (int i = 0; i < n; i++) {
			A.add(new Vector<E>());
		}

		// needed to fill the matrices
		int posEdgeId;
		double currval;

		// iterate over positive edges to build up the new cycles
		for (E e : delaunay.getPositiveEdges()) {
			// get positive edge index
			posEdgeId = adapters.get(EdgeIndex.class, e, Integer.class);
			// add edge with this index to the cycles
			for (int i = 0; i < n; i++) {
				currval = cycles.get(posEdgeId, i);
				// add the positive edge if the weight is positive else the
				// opposite edge
				E currEdge = (currval > 0) ? e : e.getOppositeEdge();
				// as many times as the weight says
				for (int count = 0; count < Math.abs(currval); count++)
					A.get(i).add(currEdge);
			}
		}

		return A;
	}
	
	/**
	 * Rotate e around v clockwise
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param e
	 * @param v
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> E getNextPrimalEdgeClockwise(E e, V v) {
		
		if (e.getStartVertex() == v) {
			return e.getPreviousEdge().getOppositeEdge();
		}
		if (e.getTargetVertex() == v) {
			return e.getOppositeEdge().getPreviousEdge();
		}
		throw new IllegalArgumentException(
				"Edge does not contain vertex in getNextEdgeClockwise()");
	}
		
	/**
	 * Rotate e around f clockwise
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param e
	 * @param f
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> E getNextDualEdgeClockwise(E e, F f) {
		
		if (e.getLeftFace() == f) {
			return e.getPreviousEdge();
		}
		if (e.getRightFace() == f) {
			return e.getOppositeEdge().getPreviousEdge().getOppositeEdge();
		}
		throw new IllegalArgumentException(
				"Edge does not contain vertex in getNextEdgeClockwise()");
	}
	
	/**
	 * There are several possibilities how an edge can be connected to the given
	 * cycle. If the surface is cut along the cycle, we can regard the cycle as
	 * two cycles actually - a left and a right one.
	 */
	private static enum EdgeStatus {
		endsAtLeftCycle,
		startsAtLeftCycle,
		endsAtRightCycle,
		startsAtRightCycle,
		liesOnLeftCycle,
		liesOnRightCycle,
		noConnection
	}

	/**
	 * Returns the edge's status, see EdgeStatus description.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param e
	 * @param edgeCycle
	 * @param vertexCycle
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> EdgeStatus getPrimalEdgeStatus(E e, List<E> edgeCycle, Set<V> vertexCycle) {
		
		if (edgeCycle.contains(e))
			return EdgeStatus.liesOnLeftCycle;
		if (edgeCycle.contains(e.getOppositeEdge()))
			return EdgeStatus.liesOnRightCycle;
		boolean outpointing = vertexCycle.contains(e.getStartVertex());
		V v = outpointing ? e.getStartVertex() : e.getTargetVertex();
		E curr = getNextPrimalEdgeClockwise(e, v);
		while (curr != e) {
			if (edgeCycle.contains(curr)) {
				if (outpointing)
					return EdgeStatus.startsAtRightCycle;
				else
					return EdgeStatus.endsAtLeftCycle;
			}
			if (edgeCycle.contains(curr.getOppositeEdge())) {
				if (outpointing)
					return EdgeStatus.startsAtLeftCycle;
				else
					return EdgeStatus.endsAtRightCycle;
			}
			curr = getNextPrimalEdgeClockwise(curr, v);
		}

		return EdgeStatus.noConnection;
	}
	
	/**
	 * Returns the dual edge's status, see EdgeStatus description.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param e
	 * @param edgeCycle
	 * @param dualVertexCycle
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> EdgeStatus getDualEdgeStatus(E e, List<E> edgeCycle, Set<F> dualVertexCycle) {
		
		if (edgeCycle.contains(e))
			return EdgeStatus.liesOnLeftCycle;
		if (edgeCycle.contains(e.getOppositeEdge()))
			return EdgeStatus.liesOnRightCycle;
		boolean outpointing = dualVertexCycle.contains(e.getRightFace());
		F f = outpointing ? e.getRightFace() : e.getLeftFace();
		E curr = getNextDualEdgeClockwise(e, f);
		while (curr != e) {
			if (edgeCycle.contains(curr)) {
				if (outpointing)
					return EdgeStatus.startsAtRightCycle;
				else
					return EdgeStatus.endsAtLeftCycle;
			}
			if (edgeCycle.contains(curr.getOppositeEdge())) {
				if (outpointing)
					return EdgeStatus.startsAtLeftCycle;
				else
					return EdgeStatus.endsAtRightCycle;
			}
			curr = getNextDualEdgeClockwise(curr, f);
		}

		return EdgeStatus.noConnection;
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

		Map<Integer, Integer> tau = getIdentificationMapOnPrimalMesh(hds, cycle);

		// Consider the surface hds being obtained by identifying two boundary
		// cycles of another surface M via the map above. On this surface we can
		// get a harmonic function with boundary conditions 0 and 1 constant on
		// the identified cycles.
		DoubleMatrix1D h = getStandardHarmonicFunctionOnPrimalMesh(hds, adapters, cycle,
				tau);

		// The differential of h can be defined on hds.
		DoubleMatrix1D dh = DoubleFactory1D.dense.make(hds.numEdges() / 2);
		Set<V> boundaryVertexSet = getPrimalVertexSet(cycle);

		// build the differences for each edge
		for (E e : hds.getPositiveEdges()) {

			int k = adapters.get(EdgeIndex.class, e, Integer.class);
			EdgeStatus status = getPrimalEdgeStatus(e, cycle, boundaryVertexSet);

			switch (status) {
			case startsAtRightCycle:
				dh.set(k,
						h.get(e.getTargetVertex().getIndex())
								- h.get(tau.get(e.getStartVertex().getIndex())));
				break;
			case endsAtRightCycle:
				dh.set(k,
						h.get(tau.get(e.getTargetVertex().getIndex()))
								- h.get(e.getStartVertex().getIndex()));
				break;
			default:
				dh.set(k,
						h.get(e.getTargetVertex().getIndex())
								- h.get(e.getStartVertex().getIndex()));
				break;
			}
		}

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

		Map<Integer, Integer> tau = getIdentificationMapOnDualMesh(hds, cycle);

		// Consider the surface hds being obtained by identifying two boundary
		// cycles of another surface M via the map above. On this surface we can
		// get a harmonic function with boundary conditions 0 and 1 constant on
		// the identified cycles.
		DoubleMatrix1D h = getStandardHarmonicFunctionOnDualMesh(hds, adapters, cycle,
				tau);

		// The differential of h can be defined on hds.
		DoubleMatrix1D dh = DoubleFactory1D.dense.make(hds.numEdges() / 2);
		Set<F> boundaryVertexSet = getDualVertexSet(cycle);

		// build the differences for each edge
		for (E e : hds.getPositiveEdges()) {

			int k = adapters.get(EdgeIndex.class, e, Integer.class);
			EdgeStatus status = getDualEdgeStatus(e, cycle, boundaryVertexSet);

			switch (status) {
			case startsAtRightCycle:
				dh.set(k,
						h.get(e.getLeftFace().getIndex())
								- h.get(tau.get(e.getRightFace().getIndex())));
				break;
			case endsAtRightCycle:
				dh.set(k,
						h.get(tau.get(e.getLeftFace().getIndex()))
								- h.get(e.getRightFace().getIndex()));
				break;
			default:
				dh.set(k,
						h.get(e.getLeftFace().getIndex())
								- h.get(e.getRightFace().getIndex()));
				break;
			}
		}

		return dh;
	}

	/**
	 * Returns a harmonic function, which is two valued on the cycle (constant 0
	 * on the left and 1 on the right side).
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param cycle
	 * @param tau
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix1D getStandardHarmonicFunctionOnPrimalMesh(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle, Map<Integer, Integer> tau) {
	
		double[] bc1 = new double[tau.size()];
		double[] bc2 = new double[tau.size()];
		for (int i = 0; i < tau.size(); i++) {
			bc2[i] = 1.;
		}
		return getPrimalHarmonicFunction(hds, adapters, cycle, tau, bc1, bc2);
		
	}

	/**
	 * Returns a harmonic function on the dual surface, which is two valued on
	 * the cycle (constant 0 on the left and 1 on the right side).
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param cycle
	 * @param tau
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix1D getStandardHarmonicFunctionOnDualMesh(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle, Map<Integer, Integer> tau) {
	
		double[] bc1 = new double[tau.size()];
		double[] bc2 = new double[tau.size()];
		for (int i = 0; i < tau.size(); i++) {
			bc2[i] = 1.;
		}
		return getDualHarmonicFunction(hds, adapters, cycle, tau, bc1, bc2);
		
	}

	/**
	 * Returns a harmonic function on the surface cut along the given cycle. The
	 * cycle splits in a left and a right side, where boundary conditions can be
	 * specified. Tau is the identification map to obtain the surface given by
	 * the hds.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param cycle
	 * @param tau
	 * @param boundaryCondition1
	 *            (left side)
	 * @param boundaryCondition2
	 *            (right side)
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix1D getPrimalHarmonicFunction(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle, Map<Integer, Integer> tau,
			double[] boundaryCondition1, double[] boundaryCondition2) {

		Set<V> vertexSet = getPrimalVertexSet(cycle);
		
		int n = hds.numVertices() + vertexSet.size();
		DoubleMatrix2D laplaceop = getPrimalLaplacian(hds, adapters, cycle, tau);

		DoubleMatrix1D x = DoubleFactory1D.dense.make(n);

		double[] bcond = new double[n];
		int bcIndex = 0;
		for (V v : vertexSet) {
			int id = v.getIndex();
			bcond[id] = boundaryCondition1[bcIndex];
			bcond[tau.get(id)] = boundaryCondition2[bcIndex];
			bcIndex++;
		}

		DoubleMatrix1D b = DoubleFactory1D.dense.make(bcond);

		// finally comment out
		// print(b,2);

		DoubleIterativeSolver solver;
		solver = new DoubleBiCG(x);
		solver = new DoubleBiCGstab(x);
		// solver= new DoubleCGLS();
		solver = new DoubleGMRES(x);

		DefaultDoubleIterationMonitor monitor = new DefaultDoubleIterationMonitor();

		// DoubleIterationMonitor monitor= new CGLSDoubleIterationMonitor();
		// configure monitor
		monitor.setMaxIterations(100000);
		monitor.setAbsoluteTolerance(eps);
		monitor.setRelativeTolerance(eps);

		monitor.setDivergenceTolerance(1);

		monitor.setNormType(Norm.Two);

		monitor.setIterationReporter(reporter);

		solver.setIterationMonitor(monitor);

		try {
			solver.solve(laplaceop, b, x);
		} catch (IterativeSolverDoubleNotConvergedException e) {
			System.err
					.println("Iterative solver failed to converge: Couldn't get harmonic function.");
			e.printStackTrace();
		}
		System.err.println();

		DoubleMatrix1D H = DoubleFactory1D.dense.make(n);

		for (int i = 0; i < n; i++) {
			H.set(i, x.get(i));
		}

		// for (Integer I : tau.keySet()) {
		// System.err.println("h+ = " + H.get(I) + " ; h- = "
		// + H.get(tau.get(I)));
		// }

		return H;
	}

	/**
	 * Returns a harmonic function on the dual surface cut along the given cycle. The
	 * cycle splits in a left and a right side, where boundary conditions can be
	 * specified. Tau is the identification map to obtain the surface given by
	 * the dual hds.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param cycle
	 * @param tau
	 * @param boundaryCondition1
	 *            (left side)
	 * @param boundaryCondition2
	 *            (right side)
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix1D getDualHarmonicFunction(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle, Map<Integer, Integer> tau,
			double[] boundaryCondition1, double[] boundaryCondition2) {

		Set<F> vertexSet = getDualVertexSet(cycle);
		int n = hds.numFaces() + vertexSet.size();
		
		DoubleMatrix2D laplaceop = getDualLaplacian(hds, adapters, cycle, tau);
		
		print(laplaceop,2);
		
		System.err.println("determinant: "+dalgebra.det(laplaceop));

		DoubleMatrix1D x = DoubleFactory1D.dense.make(n);

		double[] bcond = new double[n];
		int bcIndex = 0;
		for (F f : vertexSet) {
			int id = f.getIndex();
			bcond[id] = boundaryCondition1[bcIndex];
			bcond[tau.get(id)] = boundaryCondition2[bcIndex];
			bcIndex++;
		}

		DoubleMatrix1D b = DoubleFactory1D.dense.make(bcond);

//		 finally comment out
		 print(b,2);

		DoubleIterativeSolver solver;
		solver = new DoubleBiCG(x);
		solver = new DoubleBiCGstab(x);
		// solver= new DoubleCGLS();
		solver = new DoubleGMRES(x);

		DefaultDoubleIterationMonitor monitor = new DefaultDoubleIterationMonitor();

		// DoubleIterationMonitor monitor= new CGLSDoubleIterationMonitor();
		// configure monitor
		monitor.setMaxIterations(100000);
		monitor.setAbsoluteTolerance(eps);
		monitor.setRelativeTolerance(eps);

		monitor.setDivergenceTolerance(1);

		monitor.setNormType(Norm.Two);

		monitor.setIterationReporter(reporter);

		solver.setIterationMonitor(monitor);

		try {
			solver.solve(laplaceop, b, x);
		} catch (IterativeSolverDoubleNotConvergedException e) {
			System.err
					.println("Iterative solver failed to converge: Couldn't get harmonic function.");
			e.printStackTrace();
		}
		System.err.println();

		DoubleMatrix1D H = DoubleFactory1D.dense.make(n);

		for (int i = 0; i < n; i++) {
			H.set(i, x.get(i));
		}

		 for (Integer I : tau.keySet()) {
		 System.err.println("h+ = " + H.get(I) + " ; h- = "
		 + H.get(tau.get(I)));
		 }

		return H;
	}

	/**
	 * Returns a map which maps the indices of the vertices v0,v1,... contained
	 * in the cycle to the integers n,n+1,... , where n is the number of
	 * vertices of the surface.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param cycle
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> Map<Integer, Integer> getIdentificationMapOnPrimalMesh(
				HalfEdgeDataStructure<V, E, F> hds, List<E> cycle) {
			
		Map<Integer, Integer> map = new HashMap<Integer, Integer>(100);

		// Get the vertices contained in the cycle.
		Set<V> vertexSet = getPrimalVertexSet(cycle);
		int n = hds.numVertices();

		// Build up the map.
		for (V v:vertexSet) {
			map.put(v.getIndex(), n++);
		}
		return map;
	}

	/**
	 * Returns a map which maps the indices of the vertices v0,v1,... contained
	 * in the cycle to the integers n,n+1,... , where n is the number of
	 * vertices of the dual surface (i.e. faces).
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param cycle
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> Map<Integer, Integer> getIdentificationMapOnDualMesh(
				HalfEdgeDataStructure<V, E, F> hds, List<E> cycle) {
			
		Map<Integer, Integer> map = new HashMap<Integer, Integer>(100);

		// Get the vertices contained in the cycle.
		Set<F> vertexSet = getDualVertexSet(cycle);
		// for (F f : vertexSet) {
		// System.err.print(f.getIndex() + ", ");
		// }
		// System.err.println();
		int n = hds.numFaces();

		// Build up the map.
		for (F f:vertexSet) {
			map.put(f.getIndex(), n++);
		}
		return map;
	}
	
	/**
	 * Returns the indices of the vertices contained in the cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param cycle
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<V> getPrimalVertexSet(List<E> cycle) {
		
		Set<V> vertexSet = new HashSet<V>(cycle.size());
		for (E e:cycle) {
			vertexSet.add(e.getStartVertex());
		}
		return vertexSet;
	}
	
	/**
	 * Returns the indices of the faces contained in the dual cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param cycle
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<F> getDualVertexSet(List<E> cycle) {
		
		Set<F> vertexSet = new HashSet<F>(cycle.size());
		for (E e:cycle) {
			vertexSet.add(e.getRightFace());
		}
		return vertexSet;
	}

	/**
	 * Returns the matrix of the cotan laplace operator corresponding to the surface
	 * obtained by cutting hds along the given cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param cycle
	 * @param tau
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix2D getPrimalLaplacian(HalfEdgeDataStructure<V, E, F> hds, 
			AdapterSet adapters, List<E> cycle, Map<Integer, Integer> tau) {

		Set<V> boundaryVertexSet = getPrimalVertexSet(cycle);

		DoubleMatrix2D M= DoubleFactory2D.sparse.make(hds.numVertices()+cycle.size(),
				hds.numVertices()+cycle.size());
		
		EdgeStatus status;
		double weight;
		int i, j;

		CotanAdapter ca = new CotanAdapter();
		adapters.add(ca);
		
		for (E e : hds.getEdges()) {
			status = getPrimalEdgeStatus(e, cycle, boundaryVertexSet);
			weight = getCotanWeight(e, adapters);
			switch (status) {
			case liesOnLeftCycle:
				i = e.getStartVertex().getIndex();
				M.set(i, i, 1.);
				break;
			case liesOnRightCycle:
				i = tau.get(e.getStartVertex().getIndex());
				M.set(i, i, 1.);
				break;
			case endsAtRightCycle:
				i = e.getStartVertex().getIndex();
				j = tau.get(e.getTargetVertex().getIndex());
				M.set(i, j, weight);
				M.set(i, i, M.get(i, i) - weight);
				break;
			case endsAtLeftCycle:
				i = e.getStartVertex().getIndex();
				j = e.getTargetVertex().getIndex();
				M.set(i, i, M.get(i, i) - weight);
				M.set(i, j, weight);
				break;
			case startsAtRightCycle:
				break;
			case startsAtLeftCycle:
				break;
			case noConnection:
				i = e.getStartVertex().getIndex();
				j = e.getTargetVertex().getIndex();
				M.set(i, i, M.get(i, i) - weight);
				M.set(i, j, weight);
				break;
			}
		}
		
		adapters.remove(ca);
		return M;

	}
	
	/**
	 * Returns the matrix of the cotan laplace operator corresponding to the dual surface
	 * obtained by cutting the dual hds along the given cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param cycle
	 * @param tau
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix2D getDualLaplacian(HalfEdgeDataStructure<V, E, F> hds, 
			AdapterSet adapters, List<E> cycle, Map<Integer, Integer> tau) {
		
		Set<F> boundaryVertexSet = getDualVertexSet(cycle);

		DoubleMatrix2D M= DoubleFactory2D.sparse.make(hds.numFaces()+boundaryVertexSet.size(),
				hds.numFaces()+boundaryVertexSet.size());
		
		EdgeStatus status;
		double weight;
		int i, j;

		CotanAdapter ca = new CotanAdapter();
		adapters.add(ca);
		
		for (E e : hds.getEdges()) {
			status = getDualEdgeStatus(e, cycle, boundaryVertexSet);
			weight = 1./getCotanWeight(e, adapters);
			switch (status) {
			case liesOnLeftCycle:
				i = e.getRightFace().getIndex();
				M.set(i, i, 1.);
				// j = tau.get(e.getRightFace().getIndex());
				// M.set(j, j, 1.);
				break;
			case liesOnRightCycle:
				i = tau.get(e.getRightFace().getIndex());
				M.set(i, i, 1.);
				break;
			case endsAtRightCycle:
				i = e.getRightFace().getIndex();
				j = tau.get(e.getLeftFace().getIndex());
				M.set(i, j, weight);
				M.set(i, i, M.get(i, i) - weight);
				break;
			case endsAtLeftCycle:
				i = e.getRightFace().getIndex();
				j = e.getLeftFace().getIndex();
				M.set(i, i, M.get(i, i) - weight);
				M.set(i, j, weight);
				break;
			case startsAtRightCycle:
				break;
			case startsAtLeftCycle:
				break;
			case noConnection:
				i = e.getRightFace().getIndex();
				j = e.getLeftFace().getIndex();
				M.set(i, i, M.get(i, i) - weight);
				M.set(i, j, weight);
				break;
			}
		}
		
		adapters.remove(ca);
		return M;

	}

	private static double getCotanWeight(Edge<?, ?, ?> e, AdapterSet adapters) {
		return .5 * (adapters.get(Weight.class, e, Double.class) + adapters
				.get(Weight.class, e.getOppositeEdge(), Double.class));
	}

}
