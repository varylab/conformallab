package de.varylab.discreteconformal.util;

import java.util.List;
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
	 * Returns the matrix of the cotan laplace operator corresponding to the surface.
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
			AdapterSet adapters) {

		DoubleMatrix2D M= DoubleFactory2D.sparse.make(hds.numVertices(),
				hds.numVertices());
		
		CotanWeightAdapter ca= new CotanWeightAdapter();
		adapters.add(ca);
		
		double weight;
		int i, j;
		
		for (E e : hds.getEdges()) {
			weight = adapters.get(Weight.class, e, Double.class);
			i = e.getStartVertex().getIndex();
			j = e.getTargetVertex().getIndex();
			M.set(i, i, M.get(i, i) - weight);
			M.set(i, j, weight);
		}
		
		adapters.remove(ca);
		return M;

	}

	/**
	 * Returns the matrix of the cotan laplace operator corresponding to the dual surface.
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
			AdapterSet adapters) {

		DoubleMatrix2D M= DoubleFactory2D.sparse.make(hds.numFaces(),
				hds.numFaces());
		
		CotanWeightAdapter ca = new CotanWeightAdapter();
		adapters.add(ca);
		
		double weight;
		int i, j;
	
		for (E e : hds.getEdges()) {
			weight = 1. / adapters.get(Weight.class, e, Double.class);
			i = e.getRightFace().getIndex();
			j = e.getLeftFace().getIndex();
			M.set(i, i, M.get(i, i) - weight);
			M.set(i, j, weight);
		}
		
		adapters.remove(ca);
		return M;

	}
	
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
			AdapterSet adapters, List<E> cycle) {

		Set<V> boundaryVertexSet = EdgeUtility.getPrimalVertexSet(cycle);

		DoubleMatrix2D M= getPrimalLaplacian(hds, adapters);
		
		CotanWeightAdapter ca= new CotanWeightAdapter();
		adapters.add(ca);
		
		EdgeStatus status;
		double weight;
		int start, end;
		
		E connection= null;
		
		for (E e : hds.getEdges()) {
			status = EdgeUtility.getPrimalEdgeStatus(e, cycle, boundaryVertexSet);
			weight = adapters.get(Weight.class, e, Double.class);
			start = e.getStartVertex().getIndex();
			end = e.getStartVertex().getIndex();
			switch (status) {
			case startsAtLeftCycle:
				if(connection==null && start==cycle.get(0).getStartVertex().getIndex())
					connection= e;
				M.set(start, start, M.get(start, start) + weight);
				M.set(start, end, M.get(start, end) - weight);
				break;
			case endsAtLeftCycle:
				if(connection==null && end==cycle.get(0).getStartVertex().getIndex())
					connection= e;
				M.set(start, start, M.get(start, start) - weight);
				M.set(start, end, M.get(start, end) - weight);
				break;
			}
		}
		
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
			AdapterSet adapters, List<E> cycle) {
		
		Set<F> boundaryVertexSet = EdgeUtility.getDualVertexSet(cycle);

		DoubleMatrix2D M= getDualLaplacian(hds, adapters);
		
		CotanWeightAdapter ca= new CotanWeightAdapter();
		adapters.add(ca);
		
		EdgeStatus status;
		double weight;
		int start, end;
		
		for (E e : hds.getEdges()) {
			status = EdgeUtility.getDualEdgeStatus(e, cycle, boundaryVertexSet);
			weight = adapters.get(Weight.class, e, Double.class);
			start = e.getStartVertex().getIndex();
			end = e.getStartVertex().getIndex();
			switch (status) {
			case startsAtLeftCycle:
				M.set(start, start, M.get(start, start) + weight);
				M.set(start, end, M.get(start, end) - weight);
				break;
			case endsAtLeftCycle:
				M.set(start, start, M.get(start, start) - weight);
				M.set(start, end, M.get(start, end) - weight);
				break;
			}
		}
		
		adapters.remove(ca);
		return M;

	}

	// /**
	// * Normalizes max-norm of the rows to 1.
	// * @param M
	// */
	// private static void normalize(DoubleMatrix2D M) {
	// double curr;
	// for(int m= 0; m<M.rows(); m++){
	// double max=0;
	// for (int n = 0; n < M.columns(); n++) {
	// if(Math.abs(M.get(m, n))>Math.abs(max))
	// max=M.get(m, n);
	// }
	// if (max == 0)
	// throw new RuntimeException(
	// "ERROR: The Laplace matrix has to have a non-zero coefficient in each row!");
	// for (int n = 0; n < M.columns(); n++) {
	// curr= M.get(m, n)/max;
	// M.set(m, n, curr);
	// }
	// }
	// }

}
