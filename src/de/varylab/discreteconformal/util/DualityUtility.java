package de.varylab.discreteconformal.util;

import java.util.List;
import java.util.Set;
import java.util.Vector;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.EdgeIndex;
import de.jtem.halfedgetools.adapter.type.Weight;
import de.jtem.halfedgetools.util.HalfEdgeUtilsExtra;
import de.varylab.discreteconformal.adapter.CotanWeightAdapter;
import de.varylab.discreteconformal.util.EdgeUtility.EdgeStatus;


/**
 * Class to calculate harmonic differentials corresponding to a cycle of the
 * homology basis and holomorphic differentials in the sense of mercat.
 * 
 * By convention all the forms are saved as rows and all the cycles as columns.
 * (For private methods we always use colt matrices.)
 * 
 * @author knoeppel
 * 
 */
public class DualityUtility {

	/**
	 * Returns paths in the dual surface which are homotopic to the given
	 * ones.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param primalcycles
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<List<E>> getDualPaths(
				HalfEdgeDataStructure<V, E, F> hds, List<List<E>> primalcycles) {
		
		List<List<E>> dualcycles = new java.util.Vector<List<E>>();
		// for each cycle
		for (int i = 0; i < primalcycles.size(); i++) {
			dualcycles.add(getDualPath(hds, primalcycles.get(i)));
		}
		return dualcycles;
	}
	
	/**
	 * Returns paths in the primal surface which are homotopic to the given
	 * ones.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param dualcycles
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<List<E>> getPrimalPaths(
				HalfEdgeDataStructure<V, E, F> hds, List<List<E>> dualcycles) {
		
		List<List<E>> primalcycles = new java.util.Vector<List<E>>();
		// for each cycle
		for (int i = 0; i < dualcycles.size(); i++) {
			primalcycles.add(getPrimalPath(hds, dualcycles.get(i)));
		}
		return primalcycles;
	}
	
	/**
	 * Returns a path in the dual surface which is homotopic to the given one.
	 * 	
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param cycle
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<E> getDualPath(
				HalfEdgeDataStructure<V, E, F> hds, List<E> cycle) {
		
		List<E> dualPath = new Vector<E>();
		// get vertices contained in the primal cycle
		Set<V> vertices = EdgeUtility.getPrimalVertexSet(cycle);
		// for each vertex in the cycle
		for (V v : vertices) {
			// get the edge star
			List<E> star = HalfEdgeUtilsExtra.getEdgeStar(v);
			// and test each edge in it
			for (E e : star) {
				EdgeStatus status = EdgeUtility.getPrimalEdgeStatus(e, cycle, vertices);
				// if the edge is on the left side put it to the dual path
				// and points to the vertex
				if (status == EdgeStatus.endsAtLeftCycle)
					dualPath.add(e);
//				else if (status == EdgeStatus.startsAtLeftCycle)
//					dualPath.add(e.getOppositeEdge());
				// dualPath.add(e);
			}
		}
		EdgeUtility.removeEdgePairs(dualPath);
		return dualPath;
	}
	
	/**
	 * Returns a path in the primal surface which is homotopic to the given one.
	 * 	
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param dualcycle
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<E> getPrimalPath(
				HalfEdgeDataStructure<V, E, F> hds, List<E> dualcycle) {
		
		List<E> primalpath = new Vector<E>();
		// get vertices contained in the primal cycle
		Set<F> vertices = EdgeUtility.getDualVertexSet(dualcycle);
		// for each vertex in the cycle
		for (F f : vertices) {
			// get the edge star
			List<E> star = HalfEdgeUtilsExtra.getBoundary(f);
			// and test each edge in it
			for (E e : star) {
				EdgeStatus status = EdgeUtility.getDualEdgeStatus(e, dualcycle, vertices);
				// if the edge is on the left side put it to the dual path
				// and points to the vertex
				if (status == EdgeStatus.endsAtLeftCycle)
					primalpath.add(e);
//				else if (status == EdgeStatus.startsAtLeftCycle)
//					dualPath.add(e.getOppositeEdge());
				// dualPath.add(e);
			}
		}
		EdgeUtility.removeEdgePairs(primalpath);
		return primalpath;
	}
	
	/**
	 * Applies the Hodge star operator to an array of forms on the surface and
	 * returns the result (forms are stored as rows).
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param adapters
	 * @param forms
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix2D getDualOfPrimalForms(
				HalfEdgeDataStructure<V, E, F> delaunay, AdapterSet adapters,
				DoubleMatrix2D forms) {
		
		DoubleMatrix2D formsStar = DoubleFactory2D.dense.make(forms.rows(), forms.columns());
		DoubleMatrix1D row;
		for (int i = 0; i < formsStar.rows(); i++) {
			row= getDualOfPrimalForm(delaunay, adapters, forms.viewRow(i));
			for (int j = 0; j < forms.columns(); j++) {
				formsStar.set(i, j, row.get(j));	
			}
		}
		
		return formsStar;
	}
	
	/**
	 * Applies the Hodge star operator to an array of forms on the dual surface and
	 * returns the result (forms are stored as rows).
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param adapters
	 * @param forms
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix2D getDualOfDualForms(
				HalfEdgeDataStructure<V, E, F> delaunay, AdapterSet adapters,
				DoubleMatrix2D forms) {
		
		DoubleMatrix2D formsStar = DoubleFactory2D.dense.make(forms.rows(), forms.columns());
		DoubleMatrix1D row;
		for (int i = 0; i < formsStar.rows(); i++) {
			row= getDualOfDualForm(delaunay, adapters, forms.viewRow(i));
			for (int j = 0; j < forms.columns(); j++) {
				formsStar.set(i, j, row.get(j));	
			}
		}
		
		return formsStar;
	}
	
	/**
	 * Applies the Hodge star operator to a form on the surface and
	 * returns the result.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param adapters
	 * @param form
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix1D getDualOfPrimalForm(
				HalfEdgeDataStructure<V, E, F> delaunay, AdapterSet adapters,
				DoubleMatrix1D form) {
		
		DoubleMatrix1D formStar = DoubleFactory1D.dense.make((int)form.size());

		// for the construction of the dual forms the cotan weights are needed
		CotanWeightAdapter ca = new CotanWeightAdapter();
		adapters.add(ca);
		
		double ratio;
		int id;

		for (E e : delaunay.getPositiveEdges()) {
			ratio = adapters.get(Weight.class, e, Double.class);
			id = adapters.get(EdgeIndex.class, e, Integer.class);
			formStar.set(id, ratio * form.get(id));
		}

		adapters.remove(ca);
		
		return formStar;
	}
	
	/**
	 * Applies the Hodge star operator to a form on the dual surface and
	 * returns the result.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param adapters
	 * @param form
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix1D getDualOfDualForm(
				HalfEdgeDataStructure<V, E, F> delaunay, AdapterSet adapters,
				DoubleMatrix1D form) {
		
		DoubleMatrix1D formStar = DoubleFactory1D.dense.make((int)form.size());

		// for the construction of the dual forms the cotan weights are needed
		CotanWeightAdapter ca = new CotanWeightAdapter();
		adapters.add(ca);
		
		double ratio;
		int id;

		for (E e : delaunay.getPositiveEdges()) {
			ratio = adapters.get(Weight.class, e, Double.class);
			id = adapters.get(EdgeIndex.class, e, Integer.class);
			formStar.set(id, -form.get(id) / ratio);
		}

		adapters.remove(ca);
		
		return formStar;
	}

}
