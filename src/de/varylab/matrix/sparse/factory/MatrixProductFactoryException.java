package de.varylab.matrix.sparse.factory;

@SuppressWarnings("serial")
public class MatrixProductFactoryException extends Exception {
	private int columns, rows;
	
	@SuppressWarnings("unchecked")
	public MatrixProductFactoryException(BlockMatrixFactory A, BlockMatrixFactory B) {
		this.columns= A.numColBlocks;
		this.rows = B.numRowBlocks;
	}
	
	@Override
	public String getMessage() {
		String message = super.getMessage();
		message += "/n"
				+ "Matrix product A*B not defined (number of columns of A = "
				+ columns + " & number of rows of B = " + rows
				+ ")";
		return message;
	}
}
