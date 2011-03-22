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
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.solver.DefaultDoubleIterationMonitor;
import cern.colt.matrix.tdouble.algo.solver.DoubleBiCGstab;
import cern.colt.matrix.tdouble.algo.solver.DoubleIterationReporter;
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
			System.err.println("iteration = "+arg1+", value = "+arg0);
		}
	};
	
	private static DenseDoubleAlgebra dalgebra = new DenseDoubleAlgebra();
	private static double eps = 1E-10;
	
	/**
	 * Calculates 2*g harmonic differentials on hds. The weight Adapter is used
	 * to find a basis of the homology that are short with respect to the weight
	 * given by wa.
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
	> double[][] getHarmonicForms(
		HalfEdgeDataStructure<V, E, F> hds, 
		AdapterSet adapters,
		WeightAdapter<E> wa
	) {
		// First make clear that we are working with a delaunay triangulation.
		MappedLengthAdapter la = Delaunay.constructDelaunay(hds, adapters);

		// Get the homology basis of the surface.
		V rootV = hds.getVertex(0);
		List<List<E>> basis = getCanonicalHomologyBasis(rootV,adapters, wa);
		
		// use the private method
		return getHarmonicForms(hds, basis, adapters, la, wa);
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

		// use the private method
		return getHolomorphicForms(hds, basis, adapters, la, wa);
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
	 * Calculates a canonical homology basis.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param root
	 * @param wa
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<List<E>> getCanonicalHomologyBasis(V root, AdapterSet adapters, WeightAdapter<E> wa){
		
		// get a homology basis (not canonical)
		List<Set<E>> paths=HomologyUtility.getGeneratorPaths(root, wa);
		
		// write the sets as cycles, i.e. lists (an edge can appear more than only once)
		List<List<E>> homologyBasis= new Vector<List<E>>();
		for (Set<E> path: paths) {
			Vector<E> cycle= new Vector<E>();
			for(E e: path){
				cycle.add(e);
			}
			homologyBasis.add(cycle);
		}
		
		// dimension of the homology group
		int dimension= homologyBasis.size();

		DoubleMatrix2D S = DoubleFactory2D.sparse.make(dimension, dimension);

		// // TEST: show intersection number for initial basis
		// 
		// for (int i = 0; i < dimension; i++) {
		// for (int j = i+1; j < dimension; j++) {
		// int s = getIntersectionNumber(homologyBasis.get(i),
		// homologyBasis.get(j));
		// if (s != 0){
		// S.set(i, j, s);
		// S.set(j, i, -s);
		// }
		// }
		// }
		// print(S, 0);
		
		// number of positive edges
		int numOfPosEdges= root.getHalfEdgeDataStructure().numEdges()/2;
		int posEdgeId; double currval;
		
		// create a matrix to store the basis vectors
		DoubleMatrix2D Basis= DoubleFactory2D.dense.make(numOfPosEdges,dimension);
		
		// for each path
		for (int d= 0 ; d<dimension; d++) {
			for (E e : homologyBasis.get(d)) {
				posEdgeId = adapters.get(EdgeIndex.class, e, Integer.class);
				// get current value
				currval = Basis.get(posEdgeId,d);
				// add + or -1
				if (e.isPositive())
					currval++;
				else
					currval--;
				// set weight
				Basis.set(posEdgeId,d, currval);
			}
		}
		
		// get coordinates of a canonical basis
		DoubleMatrix2D coords = canonicalBasisCoords(new BilinearForm(S));
		
		// build up new homology basis
		
		// calculate basis vectors
		DoubleMatrix2D BPrime = dalgebra.mult(Basis, coords);

		// initialize empty list of cycles
		List<List<E>> newHomologyBasis = new Vector<List<E>>(dimension);
		for (int i = 0; i < dimension; i++) {
			newHomologyBasis.add(new Vector<E>());
		}

		// iterate over positive edges to build up the new cycles
		for (E e : root.getHalfEdgeDataStructure().getPositiveEdges()) {
			// get positive edge index
			posEdgeId = adapters.get(EdgeIndex.class, e, Integer.class);
			// add edge with this index to the cycles
			for (int i = 0; i < dimension; i++) {
				currval = BPrime.get(posEdgeId, i);
				// add the positive edge if the weight is positive else the
				// opposite edge
				E currEdge = (currval > 0) ? e : e.getOppositeEdge();
				// as many times as the weight says
				for (int count = 0; count < Math.abs(currval); count++)
					newHomologyBasis.get(i).add(currEdge);
			}
		}

		// TEST: show intersection number for canonical basis
		System.out
				.println("Intersection matrix for canonical homology basis:");
		S = DoubleFactory2D.sparse.make(dimension, dimension);
		for (int i = 0; i < dimension; i++) {
			for (int j = i + 1; j < dimension; j++) {
				int s = getIntersectionNumber(newHomologyBasis.get(i),
						newHomologyBasis.get(j));
				if (s != 0) {
					S.set(i, j, s);
					S.set(j, i, -s);
				}
			}
		}

		print(S, 0);
		
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
		int dimension= B.getDimension();

		// create a coordinate matrix
		DoubleMatrix2D coords= DoubleFactory2D.dense.make(dimension,dimension);
		for (int i = 0; i < dimension; i++) {
			coords.set(i, i, 1);
		}
		
		// initialize empty index sets
		Stack<Integer> todo= new Stack<Integer>();
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

			// // Show step result
			// System.out.println("Transformation:");
			// print(T, 2);
			// System.out.println();
			// System.out.println("Coords:");

			coords = dalgebra.mult(coords, T);

			// print(coords, 2);

//			DoubleMatrix2D A = DoubleFactory2D.dense.make(dimension, dimension);
//			for (int k = 0; k < dimension; k++) {
//				for (int l = 0; l < dimension; l++) {
//					A.set(k,
//							l,
//							B.brackets(coords.viewColumn(k),
//									coords.viewColumn(l)));
//				}
//			}

			// System.out.println();
			// System.out.println("A=");
			// print(A, 2);
			//
			// System.out.println();
		}

		DoubleMatrix2D S= DoubleFactory2D.sparse.make(dimension,dimension);
		DoubleMatrix2D finalCoords= DoubleFactory2D.dense.make(dimension,dimension);
		
		for (int k = 0; k < dimension; k++) {
			for (int l = k + 1; l < dimension; l++) {
				double val = B.brackets(coords.viewColumn(k),
						coords.viewColumn(l));
				S.set(k, l, val);
				S.set(l, k, -val);
			}
		}
		
		IntArrayList rows= new IntArrayList();
		IntArrayList cols= new IntArrayList();
		DoubleArrayList vals= new DoubleArrayList();
		S.getPositiveValues(rows, cols, vals);
		
		for(int j= 0; j<rows.size() ; j++){
			for (int i = 0; i < dimension; i++) {
				finalCoords.set(i,j ,coords.get(i, rows.get(j)));
				finalCoords.set(i,j+cols.size(),coords.get(i, cols.get(j)));
			}
		}
		
		// DoubleMatrix2D A = DoubleFactory2D.dense.make(dimension, dimension);
		// for (int k = 0; k < dimension; k++) {
		// for (int l = 0; l < dimension; l++) {
		// A.set(k,
		// l,
		// B.brackets(finalCoords.viewColumn(k),
		// finalCoords.viewColumn(l)));
		// }
		// }
		//
		// print(A,2);
		
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
		double min = 1E+16;
		double curr;
		DoubleMatrix1D col = coords.viewColumn(I);
		for (int k = 0; k < B.getDimension(); k++) {
			curr = Math.abs(B.brackets(col, coords.viewColumn(k)));
			if (curr < min && curr!=0) {
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
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<List<E>> getDualPaths(HalfEdgeDataStructure<V,E,F> hds, List<List<E>> cycles){
		List<List<E>> dualCycles= new java.util.Vector<List<E>>();
		// for each cycle
		for (int i = 0; i < cycles.size(); i++) {
			dualCycles.add(getDualPath(hds, cycles.get(i)));
		}
		return dualCycles;
	}
	
	/**
	 * Returns a path in the dual surface which is homotopic to the given one.
	 * Set	
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
	> List<E> getDualPath(HalfEdgeDataStructure<V,E,F> hds, List<E> cycle){
		List<E> dualPath = new Vector<E>();
		// get vertices contained in the primal cycle
		Set<V> vertices = getVertexSet(cycle);
		// for each vertex in the cycle
		for (V v : vertices) {
			// get the edge star
			List<E> star = HalfEdgeUtilsExtra.getEdgeStar(v);
			// and test each edge in it
			for (E e : star) {
				EdgeStatus status = getEdgeStatus(e, cycle, vertices);
				// if the edge is on the left side put it to the dual path
				if (status == EdgeStatus.endsAtLeftCycle
						|| status == EdgeStatus.startsAtLeftCycle)
					dualPath.add(e);
			}
		}
		return dualPath;
	}

	/**
	 * Returns a matrix A needed to build the normalized holomorphic
	 * differentials. The matrix has format 2g*2g. There are 2g harmonic
	 * differentials dH_j= dh_j+i*(star)dh_j. The first g rows are filled with
	 * its a-periods, the second g rows with its values on the a-cycles of dual
	 * surface. That means in the upper g rows of Ax stands the real part and
	 * in the lower the imaginary part of the a-periods of the holomorphic
	 * differential omega=sum(x_j*dH_j).
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param DH
	 * @param dualDH
	 * @param a
	 * @param dualA
	 * @param adapters
	 * @return
	 */
	private static DoubleMatrix2D getCoefficientMatrix(
			DoubleMatrix2D DH, DoubleMatrix2D dualDH, DoubleMatrix2D a,
			DoubleMatrix2D dualA) {

		// The genus equals the number of a-cycles equals the number of rows of a
		int g = a.rows();

		// get a-period-matrix of the harmonic differentials DH
		// TODO: check whether dalgbra is too slow and we should better use
		// smpdoubleblas
		DoubleMatrix2D aPeriods = dalgebra.mult(a, DH);
		DoubleMatrix2D dualAPeriods = dalgebra.mult(dualA, dualDH);

		// The matrix has to have the format 2g*2g. 
		DoubleMatrix2D M = DoubleFactory2D.dense.make(2 * g, 2 * g);

		// for each cycle and each differential
		for (int i = 0; i < g; i++) {
			for (int j = 0; j < 2 * g; j++) {
				// fill the upper g rows with a-periods of the primal mesh
				M.set(i, j, aPeriods.get(i, j));
				// and the lower g rows with the a-periods of the dual mesh
				M.set(i + g, j, dualAPeriods.get(i, j));
			}
		}

		return M;
	}

	// private static <
	// V extends Vertex<V, E, F>,
	// E extends Edge<V, E, F>,
	// F extends Face<V, E, F>
	// > DoubleMatrix2D getCoefficientMatrix(
	// double[][] dh, double[][] dhStar,
	// List<List<E>> aCycles, List<List<E>> dualACycles,
	// AdapterSet adapters){
	//
	// // The genus equals the number of a-cycles
	// int g= aCycles.size();
	//
	// // The matrix has to have the format 2g*2g. TODO: Check whether the
	// // matrix should better be defined dense.
	// DoubleMatrix2D M= DoubleFactory2D.sparse.make(2*g, 2*g);
	//
	// // for each cycle and each differential
	// for (int i = 0; i < g; i++) {
	// for (int j = 0; j < 2 * g; j++) {
	// // fill the upper g rows with a-periods of the primal mesh
	// M.set(i, j, integrateFormOverCycle(dh[j], aCycles.get(i),
	// adapters));
	// // and the lower g rows with the a-periods of the dual mesh
	// M.set(i + g, j, integrateFormOverCycle(dhStar[j], dualACycles
	// .get(i), adapters));
	// }
	// }
	//
	// return M;
	// }
	

	// /**
	// * Returns the integral of the differential over the given cycle.
	// *
	// * @param <V>
	// * @param <E>
	// * @param <F>
	// * @param form
	// * @param cycle
	// * @param adapters
	// * @return
	// */
	// private static <
	// V extends Vertex<V, E, F>,
	// E extends Edge<V, E, F>,
	// F extends Face<V, E, F>
	// > double integrateFormOverCycle(
	// double[] form,
	// List<E> cycle,
	// AdapterSet adapters){
	// // init with zero
	// double integral = 0;
	// int id;
	// // for each edge of the cycle
	// for (E e : cycle) {
	// // get global index
	// id = adapters.get(EdgeIndex.class, e, Integer.class);
	// if (e.isPositive()) // if the edge is positive add
	// integral += form[id];
	// else
	// // if negative subtract the value of the form on the edge
	// integral -= form[id];
	// }
	// return integral;
	// }

	/**
	 * Returns g cycles of the homology basis representing the a cycles. TODO:
	 * Needs a canonical basis? It's not clear whether this really works in
	 * general. Why there shall be always g such cycles for an arbitrary basis.
	 * Would be good to visualize the returned cycles on the surface.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param homologyBasis
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<List<E>> getACycles(HalfEdgeDataStructure<V,E,F> hds, List<List<E>> homologyBasis){
		
		// the genus has to be one at least
		if (homologyBasis.size() < 1)
			return null;

		// a basis contains 2g cycles
		int g = homologyBasis.size() / 2;

		// g non-intersecting cycles shall be labeled as a-cycles
		List<List<E>> aCycles = new Vector<List<E>>(g);
		aCycles.add(homologyBasis.get(0));

		// try to find g cycles, which does not intersect
		for (int i = 1; i < g; i++) {
			// try to find a new cycle with intersection number 0
			for (List<E> c : homologyBasis) {
				// if the cycle is already in the list, ignore it
				if (aCycles.contains(c))
					continue;
				// say c is non-intersecting
				boolean nonintersecting = true;
				// check the intersection number with each cycle already
				// collected
				for (List<E> a : aCycles) {
					// if the intersection number is not equal to zero, c is
					// intersecting and has to be left out
					if (getIntersectionNumber(c, a) != 0) {
						nonintersecting = false;
						break;
					}
				}
				// if c is not intersecting any of the already collected cycles
				// add it to the list
				if (nonintersecting) {
					aCycles.add(c);
					break;
				}
			}
		}

		if (aCycles.size() != g)
			throw new RuntimeException(
					"Couldn't find g non-intersecting cycles.");

		return aCycles;
	}
	
	/**
	 * Returns the intersection number of two oriented cycles.
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
	> int getIntersectionNumber(List<E> cycle1, List<E> cycle2) {
		// get the number of positive intersections of the first cycle with the
		// second
		int numPositiveIntersections = countEdgesWithStatus(cycle1, cycle2,
				EdgeStatus.startsAtLeftCycle);
		// get the number of negative intersections of the first cycle with the
		// second
		int numNegativeIntersections = countEdgesWithStatus(cycle1, cycle2,
				EdgeStatus.endsAtLeftCycle);
		return numPositiveIntersections - numNegativeIntersections;
	}

	/**
	 * Count the edges in a given list having a specified status with respect to
	 * a given cycle.
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
	> int countEdgesWithStatus(List<E> cycle, List<E> list, EdgeStatus status){
		int count = 0;
		Set<V> vertexSet= getVertexSet(cycle);
		for(E e:list){
			if(getEdgeStatus(e, cycle,vertexSet)== status)
				count++;
		}
		return count;
	}
	
	/**
	 * Applies the Hodge star operator to an array of forms on the surface and
	 * returns the result.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param adapters
	 * @param forms
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double[][] getDualForms(
		HalfEdgeDataStructure<V, E, F> delaunay,
		AdapterSet adapters,
		double[][] forms){
		
		
		double[][] formsStar = new double[forms.length][];
		for (int i = 0; i < formsStar.length; i++) {
			formsStar[i] = getDualForm(delaunay, adapters, forms[i]);
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
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double[] getDualForm(
		HalfEdgeDataStructure<V, E, F> delaunay,
		AdapterSet adapters,
		double[] form){
		
		double[] formStar = new double[form.length];

		// for the construction of the dual forms the cotan weights are needed
		CotanAdapter ca = new CotanAdapter();
		adapters.add(ca);
		
		double ratio;
		int id;

		for (E e : delaunay.getPositiveEdges()) {
			ratio = getCotanWeight(e, adapters);
			id = adapters.get(EdgeIndex.class, e, Integer.class);
			formStar[id] = ratio * form[id];
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
	> double[][] getHarmonicForms(
		HalfEdgeDataStructure<V, E, F> delaunay, 
		List<List<E>> homologyBasis,
		AdapterSet adapters, MappedLengthAdapter la,
		WeightAdapter<E> wa
	) {
		
		// For each cycle in the basis we get a corresponding harmonic
		// differential (closed but not exact, its integral differs by one on
		// the cycle).
		double[][] dh = new double[homologyBasis.size()][];

		// For each cycle construct the corresponding harmonic differential.
		adapters.add(la); // delaunay lengths
		for (int i = 0; i < homologyBasis.size(); i++) {
			dh[i] = getStandardHarmonicForm(delaunay, adapters,
					homologyBasis.get(i));
		}
		adapters.remove(la);
		
		return dh;
	}
	
	/**
	 * Returns a basis of holomorphic differentials on the surface, which is
	 * normalized with respect to the given homology basis. The surface has to
	 * be Delaunay triangulated. The forms are saved as rows. The real part of
	 * the value in the i-th row corresponds to the value of the differential on
	 * the positive oriented edge with edgeIndex i, the imaginary part to the
	 * value of its dual edge.
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
	> Complex[][] getHolomorphicForms(
		HalfEdgeDataStructure<V, E, F> delaunay,
		List<List<E>> homologyBasis,
		AdapterSet adapters,
		MappedLengthAdapter la,
		WeightAdapter<E> wa){

		// the forms are defined on the positive oriented edges
		int numPosEdges = delaunay.numEdges() / 2;
		
		// g is simply the genus of the surface
		int g = homologyBasis.size() / 2;

		// Get the harmonic differentials on the surface and its dual. The
		// format of the matrices is 2g*numEdges
		double[][] dh = getHarmonicForms(delaunay, homologyBasis, adapters, la,
				wa);
		double[][] dhStar = getDualForms(delaunay, adapters, dh);

		// We want to have matrices each column of which corresponds to one
		// harmonic differential. So we can build linear combinations of the
		// harmonic differentials simply by matrix multiplication. TODO: Check,
		// whether, transpose is perhaps too expensive?
		DoubleMatrix2D DH = dalgebra.transpose(DoubleFactory2D.dense
				.make(dh));
		DoubleMatrix2D dualDH = dalgebra
				.transpose(DoubleFactory2D.dense.make(dhStar));

		// to normalize the differentials we need the a-periods and its dual
		// cycles
		List<List<E>> acycles = getACycles(delaunay, homologyBasis);
		List<List<E>> dualACycles = getDualPaths(delaunay, acycles);

		// initialize two matrices to encode the cycles
		DoubleMatrix2D A = new SparseDoubleMatrix2D(g, numPosEdges);
		DoubleMatrix2D dualA = new SparseDoubleMatrix2D(g, numPosEdges);

		// needed to fill the matrices
		int currEdgeId;
		double currval;

		// fill the matrix of a-cycles
		for (int i = 0; i < acycles.size(); i++) {
			for (E e : acycles.get(i)) {
				currEdgeId = adapters.get(EdgeIndex.class, e, Integer.class);
				currval = A.get(i, currEdgeId);
				if (e.isPositive()) {
					A.set(i, currEdgeId, currval + 1);
				} else {
					A.set(i, currEdgeId, currval - 1);
				}
			}
		}
		
		// fill the matrix of dual a-cycles
		for (int i = 0; i < dualACycles.size(); i++) {
			for (E e : dualACycles.get(i)) {
				currEdgeId = adapters.get(EdgeIndex.class, e, Integer.class);
				currval = A.get(i, currEdgeId);
				if (e.isPositive()) {
					dualA.set(i, currEdgeId, currval + 1);
				} else {
					dualA.set(i, currEdgeId, currval - 1);
				}
			}
		}

		// an array to save the forms: The real part corresponds to the
		// differential on the primal mesh and the imaginary part to the
		// differentials on the dual mesh. 
		Complex[][] OMEGA= new Complex[g][numPosEdges];
		
		// The holomorphic forms omega are linear combinations of 2g harmonic
		// forms dh_j plus i times its duals (star)dh_j: omega= sum(x_j*(dh_j+
		// i*(star)dh_j)). They are normalized iff they satisfy a linear system.
		// A is the coefficient matrix.
		DoubleMatrix2D M = getCoefficientMatrix(DH, dualDH,A,dualA);

		// The number of harmonic forms on the surface is 2g, so we need 2g
		// coefficients.
		DoubleMatrix1D x = DoubleFactory1D.dense.make(2 * g);

		// Iterative solver
		DoubleBiCGstab solver = new DoubleBiCGstab(x);
		DefaultDoubleIterationMonitor monitor= new DefaultDoubleIterationMonitor();
		
		// configure monitor
		monitor.setMaxIterations(100000);
		monitor.setAbsoluteTolerance(eps);
		monitor.setRelativeTolerance(eps);
		monitor.setIterationReporter(reporter);
		
		solver.setIterationMonitor(monitor);

		// For each a-cycle find the holomorphic form which is 1 along this
		// cycle, i.e. gives one for the primal cycle and 0 for its dual cycle. 
		for (int i = 0; i < g; i++) {
			
			// set up the conditions 
			DoubleMatrix1D bc = DoubleFactory1D.dense.make(2 * g);
			bc.set(i, 1);
			
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
			DoubleMatrix1D omega= dalgebra.mult(DH, x);
			DoubleMatrix1D omegaStar= dalgebra.mult(dualDH, x);
			
			for (int j = 0; j < numPosEdges; j++) {
				OMEGA[i][j]= new Complex(omega.get(j),omegaStar.get(j));
			}
		}
		
		adapters.remove(la);

		return OMEGA;
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
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> E getNextEdgeClockwise(E e, V v) {
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
	> EdgeStatus getEdgeStatus(E e, List<E> edgeCycle, Set<V> vertexCycle) {
		boolean outpointing = vertexCycle.contains(e.getStartVertex());
		if (edgeCycle.contains(e))
			return EdgeStatus.liesOnLeftCycle;
		if (edgeCycle.contains(e.getOppositeEdge()))
			return EdgeStatus.liesOnRightCycle;
		V v = outpointing ? e.getStartVertex() : e.getTargetVertex();
		E curr = getNextEdgeClockwise(e, v);
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
			curr = getNextEdgeClockwise(curr, v);
		}

		return EdgeStatus.noConnection;
	}
	
	/**
	 * Returns the standard harmonic differential for the given cycle, i.e. the
	 * harmonic differential, whose integral along the cycle is 1. 
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
	> double[] getStandardHarmonicForm(HalfEdgeDataStructure<V, E, F> hds, 
		AdapterSet adapters, List<E> cycle) {

		Map<Integer, Integer> tau = getIdentificationMap(hds, cycle);

		// Consider the surface hds being obtained by identifying two boundary
		// cycles of another surface M via the map above. On this surface we can
		// get a harmonic function with boundary conditions 0 and 1 constant on
		// the identified cycles.
		double[] h = getStandardHarmonicFunction(hds, adapters, cycle, tau);

		// The differential of h can be defined on hds.
		double[] dh = new double[hds.numEdges() / 2];
		Set<V> boundaryVertexSet = getVertexSet(cycle);

		// build the differences for each edge
		for (E e : hds.getPositiveEdges()) {
			
			int k= adapters.get(EdgeIndex.class, e,Integer.class);
			EdgeStatus status= getEdgeStatus(e, cycle, boundaryVertexSet);
				
			switch (status) {
			case startsAtRightCycle:
				dh[k] = h[e.getTargetVertex().getIndex()]
						- h[tau.get(e.getStartVertex().getIndex())];
				break;
			case endsAtRightCycle:
				dh[k] = h[tau.get(e.getTargetVertex().getIndex())]
						- h[e.getStartVertex().getIndex()];
				break;
			default:
				dh[k] = h[e.getTargetVertex().getIndex()]
						- h[e.getStartVertex().getIndex()];
				break;
			}
		}

		return dh;
	}

	/**
	 * Returns a harmonic function, which is two valued on the cycle
	 * (constant 0 on the left and 1 on the right side).
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
	> double[] getStandardHarmonicFunction(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle, Map<Integer, Integer> tau) {
	
		double[] bc1 = new double[tau.size()];
		double[] bc2 = new double[tau.size()];
		for (int i = 0; i < bc2.length; i++) {
			bc2[i] = 1.;
		}
		return getHarmonicFunction(hds, adapters, cycle, tau, bc1, bc2);
		
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
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double[] getHarmonicFunction(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle, Map<Integer, Integer> tau,
			double[] boundaryCondition1, double[] boundaryCondition2) {

		int n = hds.numVertices() + cycle.size();
		DoubleMatrix2D laplaceop = getLaplacian(hds, adapters, cycle, tau);

		Set<V> vertexSet = getVertexSet(cycle);

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

		DoubleBiCGstab solver = new DoubleBiCGstab(x);
		DefaultDoubleIterationMonitor monitor= new DefaultDoubleIterationMonitor();
		
		// configure monitor
		monitor.setMaxIterations(100000);
		monitor.setAbsoluteTolerance(eps);
		monitor.setRelativeTolerance(eps);
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

		double[] H = new double[n];

		for (int i = 0; i < hds.numVertices(); i++) {
			H[i] = x.get(i);
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
	> Map<Integer, Integer> getIdentificationMap(
			HalfEdgeDataStructure<V, E, F> hds, List<E> cycle) {
		Map<Integer, Integer> map = new HashMap<Integer, Integer>(100);

		// Get the vertices contained in the cycle.
		Set<V> vertexSet = getVertexSet(cycle);
		int n = hds.numVertices();

		// Build up the map.
		int i=0;
		for (V v:vertexSet) {
			map.put(v.getIndex(), n + i);
			i++;
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
	> Set<V> getVertexSet(List<E> cycle) {
		Set<V> vertexSet = new HashSet<V>(cycle.size());
		for (E e:cycle) {
			vertexSet.add(e.getStartVertex());
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
	> DoubleMatrix2D getLaplacian(HalfEdgeDataStructure<V, E, F> hds, 
			AdapterSet adapters, List<E> cycle, Map<Integer, Integer> tau) {

		Set<V> boundaryVertexSet = getVertexSet(cycle);

		DoubleMatrix2D M= DoubleFactory2D.sparse.make(hds.numVertices()+cycle.size(),
				hds.numVertices()+cycle.size());
		
		EdgeStatus status;
		double weight;
		int i, j;

		CotanAdapter ca = new CotanAdapter();
		adapters.add(ca);
		
		for (E e : hds.getEdges()) {
			status= getEdgeStatus(e, cycle, boundaryVertexSet);
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
			case startsAtRightCycle:
				i = tau.get(e.getStartVertex().getIndex());
				j = e.getTargetVertex().getIndex();
				M.set(i, j, weight);
				M.set(i, i, M.get(i, i) - weight);
				break;
			case endsAtRightCycle:
				i = e.getStartVertex().getIndex();
				j = tau.get(e.getTargetVertex().getIndex());
				M.set(i, j, weight);
				M.set(i, i, M.get(i, i) - weight);
				break;
			case endsAtLeftCycle:
			default:
				i = e.getStartVertex().getIndex();
				j = e.getTargetVertex().getIndex();
				M.set(i, j, weight);
				M.set(i, i, M.get(i, i) - weight);
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
