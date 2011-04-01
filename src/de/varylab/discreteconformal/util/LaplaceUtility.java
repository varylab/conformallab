package de.varylab.discreteconformal.util;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Weight;
import de.varylab.discreteconformal.adapter.CotanWeightAdapter;

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
		
		for (E e : hds.getPositiveEdges()) {
			weight = adapters.get(Weight.class, e, Double.class);
			i = e.getStartVertex().getIndex();
			j = e.getTargetVertex().getIndex();
			M.set(i, i, M.get(i, i) - weight);
			M.set(j, j, M.get(j, j) - weight);
			M.set(i, j, weight);
			M.set(j, i, weight);
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
	
		for (E e : hds.getPositiveEdges()) {
			weight = 1. / adapters.get(Weight.class, e, Double.class);
			i = e.getRightFace().getIndex();
			j = e.getLeftFace().getIndex();
			M.set(i, i, M.get(i, i) - weight);
			M.set(j, j, M.get(j, j) - weight);
			M.set(i, j, weight);
			M.set(j, i, weight);
		}
		
		adapters.remove(ca);
		return M;

	}

}
