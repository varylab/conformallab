package de.varylab.discreteconformal.util;

import java.util.HashMap;
import java.util.Map;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.Weight;
import de.varylab.discreteconformal.adapter.CotanWeightAdapter;
import de.varylab.discreteconformal.adapter.MappedWeightAdapter;

/**
 * Class to calculate Laplace operator of a given half edge data structure or its dual.
 * 
 * @author knoeppel
 * 
 */
public class LaplaceUtility {
	
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> MappedWeightAdapter calculateCotanWeights(HDS hds, AdapterSet a) {
		Map<Edge<?,?,?>, Double> map = new HashMap<Edge<?,?,?>, Double>();
		for (E e : hds.getEdges()) {
			double leftcotanalpha = getLeftWeight(e, a);
			double rightcotanalpha = getLeftWeight(e.getOppositeEdge(), a);
			double w = 0.5 * (leftcotanalpha + rightcotanalpha);
			map.put(e, w);
		}
		MappedWeightAdapter result = new MappedWeightAdapter(map);
		return result;
	}
	
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double getLeftWeight(E e, AdapterSet aSet) {
		Double a = aSet.get(Length.class, e, Double.class);
		Double b = aSet.get(Length.class, e.getNextEdge(), Double.class);
		Double c = aSet.get(Length.class, e.getNextEdge().getNextEdge(), Double.class);
		if (!e.getNextEdge().getNextEdge().getNextEdge().equals(e)) {
			throw new RuntimeException("Face is not a triangle.");
		}
		double cosalpha = (b * b + c * c - a * a) / (2 * b * c);
		double cotanalpha = cosalpha / (Math.sqrt(1 - cosalpha * cosalpha));
		return cotanalpha;
	};

	/**
	 * Returns the matrix of the cotan Laplace operator corresponding to the
	 * surface.
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
		
		Double weight;
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
	 * Returns the matrix of the cotan Laplace operator corresponding to the
	 * dual surface.
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
