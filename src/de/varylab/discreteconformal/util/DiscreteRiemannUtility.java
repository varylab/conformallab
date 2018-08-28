package de.varylab.discreteconformal.util;

import static de.varylab.discreteconformal.util.CanonicalBasisUtility.getCanonicalHomologyBasis;
import static de.varylab.discreteconformal.util.DiscreteHolomorphicFormUtility.getHolomorphicFormsOnDualMesh;
import static de.varylab.discreteconformal.util.DiscreteHolomorphicFormUtility.getHolomorphicFormsOnPrimalMesh;
import static java.lang.Math.PI;

import java.util.List;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import de.jtem.blas.ComplexMatrix;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.theta.SiegelReduction;
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
	
	public static class Result {
		public Complex[][] forms;
		public ComplexMatrix periodMatrix_original;
		public ComplexMatrix periodMatrix;
		public Result(Complex[][] forms, ComplexMatrix periodMatrix_original, ComplexMatrix periodMatrix) {
			this.forms = forms;
			this.periodMatrix_original = periodMatrix_original;
			this.periodMatrix = periodMatrix;
		}
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
		WeightAdapter<E> wa
	) {
		Result r = getHolomorphicFormsAndPeriodMatrix(hds, adapters, wa);
		return r.forms;
	}
		
		public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Result getHolomorphicFormsAndPeriodMatrix(
		HalfEdgeDataStructure<V, E, F> hds, 
		AdapterSet adapters,
		WeightAdapter<E> wa
	) {		
		// Get the homology basis of the surface.
		V rootV = hds.getVertex(0);
		List<List<E>> basis = getCanonicalHomologyBasis(rootV, adapters, wa);
		DoubleMatrix2D[] omega1 = getHolomorphicFormsOnPrimalMesh(hds, basis, adapters);
		DoubleMatrix2D[] omega2 = getHolomorphicFormsOnDualMesh(hds, basis, adapters);

		// to normalize the differentials we need the a-periods and its dual
		// cycles
//		List<List<E>> acycles = CanonicalBasisUtility.getACycles(basis);
//		List<List<E>> dualacycles = DualityUtility.getDualPaths(hds,acycles);

		// write cycles to matrices
//		DoubleMatrix2D A = EdgeUtility.cyclesToMatrix(adapters, hds, acycles);
//		DoubleMatrix2D dualA = EdgeUtility.cyclesToMatrix(adapters, hds, dualacycles);

//		System.out.println("A - PERIODS:");
//		System.out.println("real part of primal construction:");
//		SimpleMatrixPrintUtility.print(dalgebra.mult(omega1[0], A), 4);
//		System.out.println();
//		System.out.println("imaginary part of primal construction:");
//		SimpleMatrixPrintUtility.print(dalgebra.mult(omega1[1], dualA), 4);
//		System.out.println();
//		System.out.println("A - PERIODS:");
//		System.out.println("real part of dual construction:");
//		SimpleMatrixPrintUtility.print(dalgebra.mult(omega2[0], dualA), 4);
//		System.out.println();
//		System.out.println("imaginary part of dual construction:");
//		SimpleMatrixPrintUtility.print(dalgebra.mult(omega2[1], A), 4);
//		System.out.println();

		// to normalize the differentials we need the a-periods and its dual
		// cycles
		List<List<E>> bcycles = CanonicalBasisUtility.getBCycles(basis);
		List<List<E>> dualbcycles = DualityUtility.getDualPaths(hds, bcycles);

		// write cycles to matrices
		DoubleMatrix2D B = EdgeUtility.cyclesToMatrix(adapters, hds, bcycles);
		DoubleMatrix2D dualB = EdgeUtility.cyclesToMatrix(adapters, hds,
				dualbcycles);

		DoubleMatrix2D primalReal = dalgebra.mult(omega1[0], B);
		DoubleMatrix2D primalImag = dalgebra.mult(omega1[1], dualB);

		DoubleMatrix2D dualReal = dalgebra.mult(omega2[0], dualB);
		DoubleMatrix2D dualImag = dalgebra.mult(omega2[1], B);

//		System.out.println("PERIOD MATRIX:");
//		System.out.println("real part of primal construction:");
//		SimpleMatrixPrintUtility.print(primalReal, 4);
//		System.out.println();
//		System.out.println("imaginary part of primal construction:");
//		SimpleMatrixPrintUtility.print(primalImag, 4);
//		System.out.println();
//		System.out.println("PERIOD MATRIX:");
//		System.out.println("real part of dual construction:");
//		SimpleMatrixPrintUtility.print(dualReal, 4);
//		System.out.println();
//		System.out.println("imaginary part of dual construction:");
//		SimpleMatrixPrintUtility.print(dualImag, 4);
//		System.out.println();

		DoubleMatrix2D realPeriods = DoubleFactory2D.dense.make(primalReal
				.rows(), primalReal.columns());
		DoubleMatrix2D imagPeriods = DoubleFactory2D.dense.make(primalImag
				.rows(), primalImag.columns());
		for (int i = 0; i < primalReal.rows(); i++) {
			for (int j = 0; j < primalReal.columns(); j++) {
				realPeriods.setQuick(i, j, .5 * (primalReal.getQuick(i, j) + dualReal
						.getQuick(i, j)));
				imagPeriods.setQuick(i, j, .5 * (primalImag.getQuick(i, j) + dualImag
						.getQuick(i, j)));
			}
		}

		DoubleMatrix2D realsymm= DoubleFactory2D.dense.make(realPeriods.rows(),realPeriods.columns());
		DoubleMatrix2D imagsymm= DoubleFactory2D.dense.make(imagPeriods.rows(),imagPeriods.columns());
		for (int i=0 ; i< realPeriods.rows();i++){
			for (int j = i; j < realPeriods.columns(); j++){
				realsymm.setQuick(i,j, realPeriods.getQuick(i, j)-realPeriods.getQuick(j, i));
				realsymm.setQuick(j,i, realPeriods.getQuick(j, i)-realPeriods.getQuick(i, j));
				imagsymm.setQuick(i,j, imagPeriods.getQuick(i, j)-imagPeriods.getQuick(j, i));
				imagsymm.setQuick(j,i, imagPeriods.getQuick(j, i)-imagPeriods.getQuick(i, j));
			}
		}
			
//		System.out.println("PERIOD MATRIX (mean):");
//		System.out.println("real part:");
//		SimpleMatrixPrintUtility.print(realPeriods, 4);
//		System.out.println();
//		System.out.println("symmetry:");
//		SimpleMatrixPrintUtility.print(realsymm);
//		System.out.println();
//		System.out.println("imaginary part:");
//		SimpleMatrixPrintUtility.print(imagPeriods, 4);
//		System.out.println();
//		System.out.println("symmetry:");
//		SimpleMatrixPrintUtility.print(imagsymm);
//		System.out.println();

		int m = omega1[0].rows();
		int n = omega1[0].columns();
		Complex[][] array = new Complex[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				array[i][j] = new Complex(omega1[0].getQuick(i, j)
						+ omega2[1].getQuick(i, j), omega1[1].getQuick(i, j)
						+ omega2[0].getQuick(i, j));
			}
		}
	
		final ComplexMatrix PeriodMatrix= new ComplexMatrix(realPeriods.toArray(),imagPeriods.toArray());
		System.out.println("Period Matrix:");
		System.out.println(PeriodMatrix);
		System.out.println("Apply Siegel Reduction:");
		SiegelReduction siegel= new SiegelReduction(PeriodMatrix);
		ComplexMatrix normalizedPeriodMatrix= siegel.getReducedPeriodMatrix();
		
		normalizedPeriodMatrix.print("Normalized Period Matrix 2PI:");

		ComplexMatrix normalizedPeriodMatrix2 = normalizedPeriodMatrix.times(new Complex(0, 1.0/(2*PI)));
		normalizedPeriodMatrix2.print("Normalized Period Matrix 1:");
		
		return new Result(array, PeriodMatrix, normalizedPeriodMatrix);
	}

}
