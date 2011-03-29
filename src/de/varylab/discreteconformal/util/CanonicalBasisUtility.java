package de.varylab.discreteconformal.util;

import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.Vector;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

/**
 * Class to calculate canonical basis of homology for a given half edge data
 * structure.
 * 
 * @author knoeppel
 * 
 */
public class CanonicalBasisUtility {

	private static DenseDoubleAlgebra dalgebra = new DenseDoubleAlgebra();
	
	/**
	 * class representing a bilinear form
	 */
	private static class BilinearForm {
		public DoubleMatrix2D S;
		public int dim;

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
				int s = EdgeUtility.getIntersectionNumberOfPrimalCycles(homologyBasis.get(i),
						homologyBasis.get(j));
				if (s != 0) {
					S.set(i, j, s);
					S.set(j, i, -s);
				}
			}
		}

		// print(S, 0);
		 
		DoubleMatrix2D Basis= EdgeUtility.cyclesToMatrix(adapters, root.getHalfEdgeDataStructure(), homologyBasis);
		
		// get coordinates of a canonical basis
		DoubleMatrix2D coords = canonicalBasisCoords(new BilinearForm(S));
		
		// build up new homology basis
		
		// calculate basis vectors
		DoubleMatrix2D BPrime = dalgebra.mult(Basis, coords);

		// convert to list representation
		List<List<E>> newHomologyBasis = EdgeUtility.matrixToCycles(adapters, root.getHalfEdgeDataStructure(), BPrime);

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
	public static DoubleMatrix2D canonicalBasisCoords(BilinearForm B) {
		
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
	public static <
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
	public static <
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

}
