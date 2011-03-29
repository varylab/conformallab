package de.varylab.discreteconformal.util;

import java.util.List;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.algorithm.triangulation.Delaunay;
import de.jtem.halfedgetools.algorithm.triangulation.MappedLengthAdapter;
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
	
	private static DenseDoubleAlgebra dalgebra = new DenseDoubleAlgebra();
	
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
		List<List<E>> basis = CanonicalBasisUtility.getCanonicalHomologyBasis(
				rootV, adapters, wa);
		List<List<E>> dualbasis = DualityUtility.getDualPaths(hds,basis);
		
		DoubleMatrix2D[] omega1= DiscreteHolomorphicFormUtility.getHolomorphicFormsOnPrimalMesh(hds, basis, adapters, la, wa);
		DoubleMatrix2D[] omega2= DiscreteHolomorphicFormUtility.getHolomorphicFormsOnDualMesh(hds, dualbasis, adapters, la, wa);
		
		// to normalize the differentials we need the a-periods and its dual
		// cycles
		List<List<E>> acycles = CanonicalBasisUtility.getACycles(basis);
		List<List<E>> dualacycles = CanonicalBasisUtility.getACycles(dualbasis);

		// write cycles to matrices
		DoubleMatrix2D A = EdgeUtility.cyclesToMatrix(adapters, hds,acycles);
		DoubleMatrix2D dualA = EdgeUtility.cyclesToMatrix(adapters, hds, dualacycles);

		System.err.println();
		System.err.println("A - PERIODS:");
		System.out.println("real part of primal construction:");
		SimpleMatrixPrintUtility.print(dalgebra.mult(omega1[0], A), 4);
		System.err.println();
		System.out.println("imaginary part of primal construction:");
		SimpleMatrixPrintUtility.print(dalgebra.mult(omega1[1], dualA), 4);
		System.err.println();
		System.err.println("A - PERIODS:");
		System.out.println("real part of dual construction:");
		SimpleMatrixPrintUtility.print(dalgebra.mult(omega2[0], dualA), 4);
		System.err.println();
		System.out.println("imaginary part of dual construction:");
		SimpleMatrixPrintUtility.print(dalgebra.mult(omega2[1], A), 4);
		System.err.println();
		
		// to normalize the differentials we need the a-periods and its dual
		// cycles
		List<List<E>> bcycles = CanonicalBasisUtility.getBCycles(basis);
		List<List<E>> dualbcycles = DualityUtility.getDualPaths(hds, bcycles);

		// write cycles to matrices
		DoubleMatrix2D B = EdgeUtility.cyclesToMatrix(adapters, hds,bcycles);
		DoubleMatrix2D dualB = EdgeUtility.cyclesToMatrix(adapters, hds, dualbcycles);

		System.err.println();
		System.err.println("PERIOD MATRIX:");
		System.out.println("real part of primal construction:");
		SimpleMatrixPrintUtility.print(dalgebra.mult(omega1[0], B), 4);
		System.err.println();
		System.out.println("imaginary part of primal construction:");
		SimpleMatrixPrintUtility.print(dalgebra.mult(omega1[1], dualB), 4);
		System.err.println();
		System.err.println("PERIOD MATRIX:");
		System.out.println("real part of dual construction:");
		SimpleMatrixPrintUtility.print(dalgebra.mult(omega2[0], dualB), 4);
		System.err.println();
		System.out.println("imaginary part of dual construction:");
		SimpleMatrixPrintUtility.print(dalgebra.mult(omega2[1], B), 4);
		System.err.println();
		
		int m = omega1[0].rows();
		int n = omega1[0].columns();
		Complex[][] array = new Complex[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				array[i][j] = new Complex(omega1[0].get(i, j)
						+ omega2[1].get(i, j), omega1[1].get(i, j)
						+ omega2[0].get(i, j));
			}
		}
		
		// use the private method
		return array; 
	}

}
