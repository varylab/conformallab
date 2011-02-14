package de.varylab.matrix.sparse.factory;

public interface Multiplier<A, B> {

	/**
	 * This method defines a multiplication of to matrix blocks. This is needed
	 * for the general product of two BlockMatrix objects.
	 * 
	 * @param a
	 * @param b
	 * @return a*b (whatever this means)
	 */
	public abstract double[][] multiplyBlocks(A a, B b);
	
	public abstract int getRowBlocksize();
	public abstract int getColBlocksize();
}
