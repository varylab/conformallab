package de.varylab.matrix.sparse.factory;

public class ScalarDoubleMatrixMultilier implements Multiplier<Double, Double>{
	
	double scalar=1.;
	
	public ScalarDoubleMatrixMultilier() {
	}
	
	public ScalarDoubleMatrixMultilier(double s){
		scalar= s;
	}
	
	@Override
	public int getColBlocksize() {
		return 1;
	}
	
	@Override
	public int getRowBlocksize() {
		return 1;
	}
	
	@Override
	public double[][] multiplyBlocks(Double a, Double b) {
		return new double[][]{{scalar*a*b}};
	}

}
