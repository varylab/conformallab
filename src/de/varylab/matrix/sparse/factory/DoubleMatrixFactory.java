package de.varylab.matrix.sparse.factory;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;

public class DoubleMatrixFactory extends BlockMatrixFactory<Double> {

	public DoubleMatrixFactory(int rows, int cols) {
		super(rows, cols, 1, 1);
	}

	public DoubleMatrixFactory(Matrix M) {
		this(M.numRows(), M.numColumns());
		setMatrix(M);
	}

	/**
	 * Sets the matrix Factory ready to produce copies of M. Note: if you want
	 * to get the copy, you have to call update first.
	 * 
	 * @param M
	 */
	public void setMatrix(Matrix M) {
		numRowBlocks = M.numRows();
		numColBlocks = M.numColumns();
		clear();
		for (MatrixEntry e : M) {
			set(e.row(), e.column(), e.get());
		}
		uptodate = false;
	}

	/**
	 * Adds the entries of the given matrix to the entries collected by the
	 * factory.
	 * 
	 * @param M
	 */
	public void addMatrix(Matrix M) {
		for (MatrixEntry e : M) {
			add(e.row(), e.column(), e.get());
		}
	}

	/**
	 * Subtracts the entries of the given matrix to the entries collected by the
	 * factory.
	 * 
	 * @param M
	 */
	public void subtractMatrix(Matrix M) {
		for (MatrixEntry e : M) {
			add(e.row(), e.column(), -e.get());
		}
		uptodate = false;
	}

	/**
	 * Sets the matrix factory ready to produce copies of M transposed. Note: if
	 * you want to get the copy, you have to call update first.
	 * 
	 * @param M
	 */
	public void setMatrixTransposed(Matrix M) {
		numRowBlocks = M.numColumns();
		numColBlocks = M.numRows();
		clear();
		for (MatrixEntry e : M) {
			set(e.column(), e.row(), e.get());
		}
		uptodate = false;
	}

	@Override
	protected Double add(Double current, Double toAdd) {
		return current + toAdd;
	}

	@Override
	protected double[][] getRealMatrixRepresentation(Double entry) {
		return new double[][] { { entry } };
	}
}
