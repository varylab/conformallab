package de.varylab.matrix.sparse.factory;

import no.uib.cipr.matrix.sparse.CompRowMatrix;

public class BlockMatrix extends CompRowMatrix {

	private int numRowBlocks = 0, numColBlocks = 0;
	private int rowBlocksize = 0, colBlocksize = 0;

	public BlockMatrix(int numRowBlocks, int numColBlocks, int rowBlocksize,
			int colBlocksize, int[][] nonZeroBlockStructure) {
		super(numRowBlocks*rowBlocksize, numColBlocks*colBlocksize,
				getRealNonZeroStructure(numRowBlocks, numColBlocks,
						rowBlocksize, colBlocksize, nonZeroBlockStructure));
		this.rowBlocksize = rowBlocksize;
		this.colBlocksize = colBlocksize;
		this.numRowBlocks = numRowBlocks;
		this.numColBlocks = numColBlocks;
	}

	private static int[][] getRealNonZeroStructure(int nrb, int ncb, int rbs,
			int cbs, int[][] nz) {
		int[][] rnz = new int[nrb * rbs][];
		for (int i = 0; i < nrb; i++) {
			for (int j = 0; j < rbs; j++) {
				rnz[rbs * i + j] = new int[nz[i].length * cbs];
				for (int k = 0; k < nz[i].length; k++) {
					for (int l = 0; l < cbs; l++) {
						rnz[rbs * i + j][cbs * k + l] = cbs * nz[i][k] + l;
					}
				}
			}
		}
		return rnz;
	}
	
	public double[][] getBlock(int i, int j){
		double[][] m=new double[rowBlocksize][colBlocksize];
		for (int k = 0; k < rowBlocksize; k++) {
			for (int l = 0; l < colBlocksize; l++) {
				m[k][l]=get(rowBlocksize*i+k,colBlocksize*j+l);
			}
		}
		return m;
	}

	public int getNumRowBlocks() {
		return numRowBlocks;
	}

	public int getNumColBlocks() {
		return numColBlocks;
	}

	public int getRowBlocksize() {
		return rowBlocksize;
	}

	public int getColBlocksize() {
		return colBlocksize;
	}

}
