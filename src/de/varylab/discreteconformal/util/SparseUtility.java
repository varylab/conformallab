package de.varylab.discreteconformal.util;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.functional.Functional;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.PETSc;


public class SparseUtility {
	
	public static int[] getPETScNonZeros(int[][] sparseStucture){
		int [] nnz = new int[sparseStucture.length];
		for(int i = 0; i < nnz.length; i++){
			nnz[i] = sparseStucture[i].length;
		}
		return nnz;
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Mat getHessianTemplate(Functional<V, E, F> f, HDS hds){
		int dim = f.getDimension(hds);
		int[][] sparceStructure = f.getNonZeroPattern(hds);
		int[] nonZeros = SparseUtility.getPETScNonZeros(sparceStructure);
		Mat H = Mat.createSeqAIJ(dim, dim, PETSc.PETSC_DEFAULT, nonZeros);
		H.assemble();
		return H;
	}
	
}
