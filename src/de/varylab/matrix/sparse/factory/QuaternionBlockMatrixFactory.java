package de.varylab.matrix.sparse.factory;

import de.jreality.math.Quaternion;

public class QuaternionBlockMatrixFactory extends
		BlockMatrixFactory<Quaternion> {

	public QuaternionBlockMatrixFactory(int numRows, int numCols) {
		super(numRows, numCols, 4, 4);
	}

	@Override
	protected Quaternion add(Quaternion target, Quaternion toAdd) {
		double[] v = target.asDouble();
		return new Quaternion(v[0] + toAdd.re, v[1] + toAdd.x, v[2] + toAdd.y,
				v[3] + toAdd.z);
	}

	@Override
	protected double[][] getRealMatrixRepresentation(Quaternion entry) {
		return MatrixRepresentationFactory.makeDoubleArrayRepresentation(entry);
	}

	public static void main(String[] args) {
		System.out.println("TEST METHODS:");

		System.err.println("Set and add for the first time:");
		System.err.println();
		QuaternionBlockMatrixFactory qmf = new QuaternionBlockMatrixFactory(10,
				4);
		qmf.set(0, 0, new Quaternion(1, 2, 3, 4));
		qmf.set(3, 3, new Quaternion(1, 2, 3, 4));
		qmf.add(3, 3, new Quaternion(10, 2, 3, 4));
		qmf.add(9, 2, new Quaternion(10, 2, 3, 4));
		qmf.set(0, 1, new Quaternion());
		qmf.update();
		PrintUtils.print(qmf.getMatrix());
		System.out.println(qmf.getMatrix());
		System.out.println(QuaternionMatrixUtils
				.isQuaternionic(qmf.getMatrix()));
		System.err.println("Set and add for the second time:");
		System.err.println();
		qmf.set(0, 0, new Quaternion(-1, -8, -3, -4));
		qmf.set(3, 3, new Quaternion(-1, -8, -3, -4));
		qmf.add(3, 3, new Quaternion(-10, -8, -3, -4));
		qmf.add(9, 2, new Quaternion(-10, -8, -3, -4));
		qmf.set(0, 1, new Quaternion());
		qmf.update();
		PrintUtils.print(qmf.getMatrix());
		System.out.println(qmf.getMatrix());
		System.out.println(QuaternionMatrixUtils
				.isQuaternionic(qmf.getMatrix()));
		System.err.println("Call clear and set and add:");
		System.err.println();
		qmf.clear();
		qmf.add(3, 3, new Quaternion(-10, -8, -3, -4));
		qmf.add(9, 2, new Quaternion(-10, -8, -3, -4));
		qmf.set(0, 1, new Quaternion());
		qmf.update();
		PrintUtils.print(qmf.getMatrix());
		System.out.println(qmf.getMatrix());
		System.out.println(QuaternionMatrixUtils
				.isQuaternionic(qmf.getMatrix()));
		System.err
				.println("Change matrix format (bigger) and set entry in the matrix of new format:");
		System.err.println();
		qmf.setFormat(20, 5);
		qmf.set(19, 3, new Quaternion(-10, -8, -3, -4));
		qmf.update();
		PrintUtils.print(qmf.getMatrix());
		System.out.println(qmf.getMatrix());
		System.out.println(QuaternionMatrixUtils
				.isQuaternionic(qmf.getMatrix()));
		System.err.println("Change the matrix format (smaller):");
		System.err.println();
		qmf.setFormat(2, 20);
		qmf.update();
		PrintUtils.print(qmf.getMatrix());
		System.out.println(qmf.getMatrix());
		System.out.println(QuaternionMatrixUtils
				.isQuaternionic(qmf.getMatrix()));

		qmf.clear();
		qmf.setFormat(1, 1);
		qmf.add(0, 0, new Quaternion(-10, -8, -3, -4));
		qmf.add(0, 0, new Quaternion(-10, -8, -3, -4));
		qmf.update();
		PrintUtils.print(qmf.getMatrix());
		System.out.println(qmf.getMatrix());
		System.out.println(QuaternionMatrixUtils
				.isQuaternionic(qmf.getMatrix()));

		qmf.set(0, 0, new Quaternion());
		qmf.update();
		PrintUtils.print(qmf.getMatrix());
		System.out.println(qmf.getMatrix());
		System.out.println(QuaternionMatrixUtils
				.isQuaternionic(qmf.getMatrix()));

	}

}
