package de.varylab.matrix.sparse.factory;

import java.util.Iterator;
import java.util.Map;

import de.jreality.math.Quaternion;

public class BlockMatrixProductFactory extends BlockMatrixFactory<double[][]> {

	@SuppressWarnings("unchecked")
	private Multiplier mult;
	@SuppressWarnings("unchecked")
	private BlockMatrixFactory A = null, B = null;

	@SuppressWarnings("unchecked")
	public BlockMatrixProductFactory(Multiplier multilier) {
		super();
		mult = multilier;
	}

	@SuppressWarnings("unchecked")
	public BlockMatrixProductFactory(Multiplier multilier,
			BlockMatrixFactory A, BlockMatrixFactory B) {
		super();
		mult = multilier;
		setMatrixFactories(A, B);

	}

	/**
	 * Checks whether two matrices A, B can be multiplied, if this is possible
	 * it sets the matrices. Unless the update method is called this does
	 * nothing.
	 * 
	 * @param A
	 * @param B
	 * @throws MatrixProductFactoryException
	 */
	@SuppressWarnings("unchecked")
	public void setMatrixFactories(BlockMatrixFactory A, BlockMatrixFactory B) {
		check(A, B);
		this.A = A;
		this.B = B;
	}

	@SuppressWarnings("unchecked")
	private void check(BlockMatrixFactory A, BlockMatrixFactory B) {
		A.update();
		B.update();
		if (A.getNumColBlocks() != B.getNumRowBlocks())
			try {
				throw new MatrixProductFactoryException(A, B);
			} catch (MatrixProductFactoryException e) {
				e.printStackTrace();
			}
		numRowBlocks = A.getNumRowBlocks();
		numColBlocks = B.getNumColBlocks();

		rowBlocksize = mult.getRowBlocksize();
		colBlocksize = mult.getColBlocksize();

		nonZeroRowBlockStructureUptodate = false;
		uptodate = false;
	}

	public void update() {
		if (uptodate)
			return;
		updateProductEntries();
		super.update();
	}

	@SuppressWarnings("unchecked")
	public void updateProductEntries() {
		check(A, B);
		nonZeroRowBlockLength = new int[numRowBlocks];
		Iterator iterA = A.matrixEntries.entrySet().iterator();
		while (iterA.hasNext()) {
			Map.Entry a = (Map.Entry) iterA.next();
			Iterator iterB = B.matrixEntries.entrySet().iterator();
			MatrixIndex ida = (MatrixIndex) a.getKey();
			while (iterB.hasNext()) {
				Map.Entry b = (Map.Entry) iterB.next();
				MatrixIndex idb = (MatrixIndex) b.getKey();
				if (ida.col == idb.row)
					add(ida.row, idb.col, mult.multiplyBlocks(a.getValue(), b
							.getValue()));
			}

		}
	}

	/**
	 * Returns the number of the rows of the product matrix.
	 * 
	 * @return #rows
	 */
	public int getNumRowBlocks() {
		return numRowBlocks;
	}

	/**
	 * Returns the number of columns of the product matrix.
	 * 
	 * @return #columns
	 */
	public int getNumColBlocks() {
		return numColBlocks;
	}

	@Override
	protected double[][] add(double[][] current, double[][] toAdd) {
		double[][] sum = new double[rowBlocksize][colBlocksize];
		for (int i = 0; i < rowBlocksize; i++)
			for (int j = 0; j < colBlocksize; j++)
				sum[i][j] = current[i][j] + toAdd[i][j];
		return sum;
	}

	@Override
	protected double[][] getRealMatrixRepresentation(double[][] entry) {
		return entry;
	}

	/**
	 * main method to test some of the methods implemented.
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("Test the product factory!!!");
		BlockMatrixProductFactory mpf = new BlockMatrixProductFactory(
				new Multiplier<Quaternion, Double>() {
					@Override
					public double[][] multiplyBlocks(Quaternion a, Double b) {
						Quaternion Q = new Quaternion(a.re * b, a.x * b, a.y
								* b, a.z * b);
						return MatrixRepresentationFactory
								.makeDoubleArrayRepresentation(Q);
					}

					@Override
					public int getColBlocksize() {
						return 4;
					}

					@Override
					public int getRowBlocksize() {
						return 4;
					}
				});
		QuaternionBlockMatrixFactory A = new QuaternionBlockMatrixFactory(5, 4);
		A.add(0, 0, new Quaternion());
		A.add(0, 3, new Quaternion());
		DoubleMatrixFactory B = new DoubleMatrixFactory(4, 3);
		B.set(0, 0, 10.);
		B.set(0, 2, 10.);
		B.set(3, 1, 40.);

		mpf.setMatrixFactories(A, B);
		mpf.update();

		System.out.println("A=");
		PrintUtils.print(A.getMatrix());
		System.out.println();
		System.out.println("B=");
		PrintUtils.print(B.getMatrix());
		System.out.println();
		System.out.println("A*B=");
		PrintUtils.print(mpf.getMatrix());
		System.out.println();
	}

}
