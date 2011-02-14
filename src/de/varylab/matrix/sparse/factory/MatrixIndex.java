package de.varylab.matrix.sparse.factory;

@SuppressWarnings("unchecked")
public class MatrixIndex implements Comparable{
	
	public int row, col;

	/**
	 * Creates MatrixIndex with the specified indices.
	 * 
	 * @param row
	 * @param col
	 */
	public MatrixIndex(int row, int col) {
		this.row = row;
		this.col = col;
	}

	/**
	 * If the object is an instance of MatrixIndex: If (this.row<arg0.row)
	 * or (this.row==arg0.row and this.col<arg0.col) the method returns -1.
	 * If (this.row>arg0.row) or (this.row==arg0.row and this.col>arg0.col)
	 * the method returns 1. If (this.row==arg0.row and this.col==arg0.col)
	 * the method returns 0.
	 * 
	 * If the Object is not an instance of MatrixIndex it return -2;
	 */
	@Override
	public int compareTo(Object arg0) {
		if (arg0 instanceof MatrixIndex) {
			if (this.row == ((MatrixIndex) arg0).row) {
				if (this.col == ((MatrixIndex) arg0).col)
					return 0;
				if (this.col < ((MatrixIndex) arg0).col)
					return -1;
				if (this.col > ((MatrixIndex) arg0).col)
					return 1;
			}
			if (this.row < ((MatrixIndex) arg0).row)
				return -1;
			if (this.row > ((MatrixIndex) arg0).row)
				return 1;
		}
		return -2;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof MatrixIndex))
			return false;
		return 0 == compareTo(obj);
	}

	@Override
	public int hashCode() {
		return (row * 7) ^ col;
	}

}
