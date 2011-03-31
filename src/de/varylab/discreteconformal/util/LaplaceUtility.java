package de.varylab.discreteconformal.util;

import java.util.List;
import java.util.Map;
import java.util.Set;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Weight;
import de.varylab.discreteconformal.adapter.CotanWeightAdapter;
import de.varylab.discreteconformal.util.EdgeUtility.EdgeStatus;


/**
 * Class to calculate Laplace operator of a given half edge data structure or its dual.
 * 
 * @author knoeppel
 * 
 */
public class LaplaceUtility {
	
	/**
	 * Returns the matrix of the cotan laplace operator corresponding to the surface
	 * obtained by cutting hds along the given cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param cycle
	 * @param tau
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix2D getPrimalLaplacian(HalfEdgeDataStructure<V, E, F> hds, 
			AdapterSet adapters, List<E> cycle, Map<Integer, Integer> tau) {

		Set<V> boundaryVertexSet = EdgeUtility.getPrimalVertexSet(cycle);

		DoubleMatrix2D M= DoubleFactory2D.sparse.make(hds.numVertices()+tau.size(),
				hds.numVertices()+tau.size());
		
		CotanWeightAdapter ca= new CotanWeightAdapter();
		adapters.add(ca);
		
		EdgeStatus status;
		double weight;
		int i, j;
		
		for (E e : hds.getEdges()) {
			status = EdgeUtility.getPrimalEdgeStatus(e, cycle, boundaryVertexSet);
			weight = adapters.get(Weight.class, e, Double.class);
			switch (status) {
			case liesOnLeftCycle:
				i = e.getStartVertex().getIndex();
				M.set(i, i, 1.);
				break;
			case liesOnRightCycle:
				i = tau.get(e.getStartVertex().getIndex());
				M.set(i, i, 1.);
				break;
			case endsAtRightCycle:
				i = e.getStartVertex().getIndex();
				j = tau.get(e.getTargetVertex().getIndex());
				M.set(i, j, weight);
				M.set(i, i, M.get(i, i) - weight);
				break;
			case endsAtLeftCycle:
				i = e.getStartVertex().getIndex();
				j = e.getTargetVertex().getIndex();
				M.set(i, i, M.get(i, i) - weight);
				M.set(i, j, weight);
				break;
			case startsAtRightCycle:
				break;
			case startsAtLeftCycle:
				break;
			case noConnection:
				i = e.getStartVertex().getIndex();
				j = e.getTargetVertex().getIndex();
				M.set(i, i, M.get(i, i) - weight);
				M.set(i, j, weight);
				break;
			}
		}
		
		normalize(M);
		
		adapters.remove(ca);
		return M;

	}
	
	/**
	 * Returns the matrix of the cotan laplace operator corresponding to the dual surface
	 * obtained by cutting the dual hds along the given cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param cycle
	 * @param tau
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix2D getDualLaplacian(HalfEdgeDataStructure<V, E, F> hds, 
			AdapterSet adapters, List<E> cycle, Map<Integer, Integer> tau) {
		
		Set<F> boundaryVertexSet = EdgeUtility.getDualVertexSet(cycle);

		DoubleMatrix2D M= DoubleFactory2D.sparse.make(hds.numFaces()+tau.size(),
				hds.numFaces()+tau.size());
		
		CotanWeightAdapter ca = new CotanWeightAdapter();
		adapters.add(ca);
		
		EdgeStatus status;
		double weight;
		int i, j;
	
		for (E e : hds.getEdges()) {
			status = EdgeUtility.getDualEdgeStatus(e, cycle, boundaryVertexSet);
			weight = 1./adapters.get(Weight.class, e, Double.class);
			switch (status) {
			case liesOnLeftCycle:
				i = e.getRightFace().getIndex();
				M.set(i, i, 1.);
				break;
			case liesOnRightCycle:
				// TODO: still not fixed completely: there are still some null
				// pointer exceptions (i guess, if the mesh is not fine enough,
				// see genus 2 coarse)
				j = tau.get(e.getRightFace().getIndex());
				M.set(j, j, 1.);
				break;
			case endsAtRightCycle:
				i = e.getRightFace().getIndex();
				j = tau.get(e.getLeftFace().getIndex());
				M.set(i, j, weight);
				M.set(i, i, M.get(i, i) - weight);
				break;
			case endsAtLeftCycle:
				i = e.getRightFace().getIndex();
				j = e.getLeftFace().getIndex();
				M.set(i, i, M.get(i, i) - weight);
				M.set(i, j, weight);
				break;
			case startsAtRightCycle:
				break;
			case startsAtLeftCycle:
				break;
			case noConnection:
				i = e.getRightFace().getIndex();
				j = e.getLeftFace().getIndex();
				M.set(i, i, M.get(i, i) - weight);
				M.set(i, j, weight);
				break;
			}
		}
		
		normalize(M);
		
		adapters.remove(ca);
		return M;

	}

	/**
	 * Normalizes max-norm of the rows to 1.
	 * @param M
	 */
	private static void normalize(DoubleMatrix2D M) {
		double curr;
		for(int m= 0; m<M.rows(); m++){
			double max=0;
			for (int n = 0; n < M.columns(); n++) {
				if(Math.abs(M.get(m, n))>Math.abs(max))
					max=M.get(m, n);
			}
			if (max == 0)
				throw new RuntimeException(
						"ERROR: The Laplace matrix has to have a non-zero coefficient in each row!");
			for (int n = 0; n < M.columns(); n++) {
				curr= M.get(m, n)/max;
				M.set(m, n, curr);
			}
		}
	}

}
