package de.varylab.discreteconformal.util;


public class SparseUtility {
	
	public static int[] getPETScNonZeros(int[][] sparseStucture){
		int [] nnz = new int[sparseStucture.length];
		for(int i = 0; i < nnz.length; i++){
			nnz[i] = sparseStucture[i].length;
		}
		return nnz;
	}
	
}
