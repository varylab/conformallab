package de.varylab.matrix.sparse.factory;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import no.uib.cipr.matrix.Matrix;

public abstract class BlockMatrixFactory<E> {

	// class variables ---------------------------------------------------------

	protected Map<MatrixIndex, E> matrixEntries = new HashMap<MatrixIndex, E>();

	protected int rowBlocksize, colBlocksize;
	protected int[] nonZeroRowBlockLength;

	protected int numRowBlocks = 0, numColBlocks = 0;
	protected int[][] nonZeroRowBlockStructure;

	protected boolean nonZeroRowBlockStructureUptodate = false;
	protected BlockMatrix matrix = null;

	protected boolean uptodate = false;

	// constructor -------------------------------------------------------------

	protected BlockMatrixFactory() {
	}

	public BlockMatrixFactory(int numRowBlocks, int numColBlocks,
			int rowBlocksize, int colBlocksize) {
		setFormat(numRowBlocks, numColBlocks, rowBlocksize, colBlocksize);
	}

	// reorganizing structure methods ------------------------------------------

	/**
	 * Removes all entries of the matrix.
	 */
	public void clear() {
		matrixEntries.clear();
		nonZeroRowBlockLength = new int[numRowBlocks];
		nonZeroRowBlockStructureUptodate = false;
		uptodate = false;
	}

	/**
	 * Is called when the matrix format was changed.
	 */
	protected void cleanup() {
		int[] newNonZeroRowBlockLength = new int[numRowBlocks];
		if (nonZeroRowBlockLength != null)
			for (int i = 0; i < Math.min(newNonZeroRowBlockLength.length,
					nonZeroRowBlockLength.length); i++) {
				newNonZeroRowBlockLength[i] = nonZeroRowBlockLength[i];
			}
		nonZeroRowBlockLength = newNonZeroRowBlockLength;
		nonZeroRowBlockStructureUptodate = false;

		encompassEntries();
		uptodate = false;
	}

	/**
	 * Is called when the matrix format was changed. Removes all entries which
	 * are outside the matrix bounds after resize.
	 */
	@SuppressWarnings("unchecked")
	private void encompassEntries() {
		List<MatrixIndex> ids = new Vector<MatrixIndex>(matrixEntries.size());
		for (Map.Entry e : matrixEntries.entrySet()) {
			MatrixIndex id = (MatrixIndex) e.getKey();
			if (id.col >= numColBlocks || id.row >= numRowBlocks) {
				ids.add(id);
			}
		}
		MatrixIndex currid;
		for (int i = 0; i < ids.size(); i++) {
			currid = ids.get(i);
			matrixEntries.remove(currid);
			if (currid.col >= numColBlocks && currid.row < numRowBlocks)
				nonZeroRowBlockLength[currid.row]--;
		}
		nonZeroRowBlockStructureUptodate = false;
		uptodate = false;
	}

	// update methods ---------------------------------------------------------

	/**
	 * Updates the CompRowMatrix. Has to be called before the changes are
	 * applied to the CompRowMatrix. Take care: Initializes a new sparse matrix
	 * from the gives entries, i.e. the pointer changes.
	 */
	@SuppressWarnings("unchecked")
	public void update() {
		if (uptodate)
			return;
		if (!nonZeroRowBlockStructureUptodate)
			updateNonZeroStructure();
		matrix = new BlockMatrix(numRowBlocks, numColBlocks, rowBlocksize,
				colBlocksize, nonZeroRowBlockStructure);
		MatrixIndex id;
		E entry;
		for (Map.Entry e : matrixEntries.entrySet()) {
			id = (MatrixIndex) e.getKey();
			entry = (E) e.getValue();
			writeIntoRealMatrix(id.row, id.col, entry);
		}
		uptodate = true;
	}

	/**
	 * Updates the non zero block structure needed for building a block matrix.
	 */
	@SuppressWarnings("unchecked")
	private void updateNonZeroStructure() {
		nonZeroRowBlockStructure = new int[numRowBlocks][];
		for (int i = 0; i < numRowBlocks; i++) {
			nonZeroRowBlockStructure[i] = new int[nonZeroRowBlockLength[i]];
		}
		MatrixIndex id;
		int[] numCurrEntryOfRow = new int[numRowBlocks];
		for (Map.Entry e : matrixEntries.entrySet()) {
			id = (MatrixIndex) e.getKey();
			nonZeroRowBlockStructure[id.row][numCurrEntryOfRow[id.row]] = id.col;
			numCurrEntryOfRow[id.row]++;
		}
		nonZeroRowBlockStructureUptodate = false;
	}

	// write real matrix method ------------------------------------------------

	protected void writeIntoRealMatrix(int i, int j, E entry) {
		double[][] block = getRealMatrixRepresentation(entry);
		for (int k = 0; k < rowBlocksize; k++)
			for (int l = 0; l < colBlocksize; l++)
				matrix.set(rowBlocksize * i + k, colBlocksize * j + l,
						block[k][l]);
	}

	// getters and setters -----------------------------------------------------

	/**
	 * Returns the quaternion valued matrix. In fact, this method returns a CRS
	 * matrix. Note: That your changes are applied to the matrix, you have to
	 * call the update method.
	 * 
	 * @return matrix
	 */
	public Matrix getMatrix() {
		return matrix;
	}

	/**
	 * Returns the number of the rows of the quaternion valued matrix.
	 * 
	 * @return number of rows
	 */
	public int getNumRowBlocks() {
		return numRowBlocks;
	}

	/**
	 * Returns the number of the columns of the quaternion valued matrix.
	 * 
	 * @return number of columns
	 */
	public int getNumColBlocks() {
		return numColBlocks;
	}

	/**
	 * Returns the number of the rows of the real matrix.
	 * 
	 * @return number of real rows
	 */
	public int getRealNumRows() {
		return numRowBlocks * rowBlocksize;
	}

	/**
	 * Returns the number of the columns of the real matrix.
	 * 
	 * @return number of real columns
	 */
	public int getRealNumCols() {
		return numColBlocks * colBlocksize;
	}

	/**
	 * Sets the format of the matrix which shall be produced.
	 * 
	 * @param numRowBlocks
	 * @param numColBlocks
	 * @param rowBlocksize
	 * @param colBlocksize
	 */
	protected void setFormat(int numRowBlocks, int numColBlocks,
			int rowBlocksize, int colBlocksize) {
		this.numRowBlocks = numRowBlocks;
		this.numColBlocks = numColBlocks;
		this.rowBlocksize = rowBlocksize;
		this.colBlocksize = colBlocksize;
		cleanup();
	}

	public void setFormat(int numRowBlocks, int numColBlocks) {
		setFormat(numRowBlocks, numColBlocks, rowBlocksize, colBlocksize);
	}

	public void setBlockFormat(int rowBlocksize, int colBlocksize) {
		setFormat(numRowBlocks, numColBlocks, rowBlocksize, colBlocksize);
	}

	// set, add, remove --------------------------------------------------------

	public void set(int i, int j, E entry) {
		MatrixIndex id = new MatrixIndex(i, j);
		if (!matrixEntries.containsKey(id)) {
			nonZeroRowBlockLength[i]++;
			nonZeroRowBlockStructureUptodate = false;
		}
		matrixEntries.put(id, entry);
		uptodate = false;
	}

	public void add(int i, int j, E entry) {
		MatrixIndex id = new MatrixIndex(i, j);
		E oldEntry = matrixEntries.get(id);
		if (oldEntry == null) {
			nonZeroRowBlockLength[i]++;
			nonZeroRowBlockStructureUptodate = false;
			matrixEntries.put(id, entry);
		} else {
			matrixEntries.put(id, add(oldEntry, entry));
		}
		uptodate = false;
	}

	public void remove(int i, int j) {
		MatrixIndex id = new MatrixIndex(i, j);
		matrixEntries.remove(id);
		if (nonZeroRowBlockLength[i] > 0)
			nonZeroRowBlockLength[i]--;
		else
			System.err.println("Error: all row entries are zero already.");
		nonZeroRowBlockStructureUptodate = false;
		uptodate = false;
	}

	// abstract methods --------------------------------------------------------

	/**
	 * Has to be implemented for the writing the entries of the hash map to the
	 * real (CRS) matrix. The the array sure has to have the correct dimension.
	 * That means here [rowBlocksize][colBlocksize].
	 * 
	 * @param entry
	 * @return the corresponding matrix block
	 */
	abstract protected double[][] getRealMatrixRepresentation(E entry);

	/**
	 * Has to be implemented for the use in the public add method. Shall add the
	 * toAdd to the target.
	 * 
	 * @param current
	 * @param toAdd
	 */
	abstract protected E add(E current, E toAdd);

}
