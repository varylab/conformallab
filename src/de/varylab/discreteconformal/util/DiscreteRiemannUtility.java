package de.varylab.discreteconformal.util;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import no.uib.cipr.matrix.sparse.GMRES;
import no.uib.cipr.matrix.sparse.IterativeSolver;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import no.uib.cipr.matrix.sparse.OutputIterationReporter;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.EdgeIndex;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.Weight;
import de.jtem.halfedgetools.adapter.type.generic.Position3d;
import de.jtem.halfedgetools.algorithm.triangulation.Delaunay;
import de.jtem.halfedgetools.algorithm.triangulation.DelaunayLengthAdapter;
import de.jtem.halfedgetools.util.HalfEdgeUtilsExtra;
import de.varylab.discreteconformal.util.Search.WeightAdapter;
import de.varylab.matrix.sparse.factory.BlockMatrixFactory;
import de.varylab.matrix.sparse.factory.DoubleMatrixFactory;

public class DiscreteRiemannUtility {

	@Weight
	private static class CotanAdapter extends AbstractAdapter<Double> {

		public CotanAdapter() {
			super(Double.class, true, false);
		}

		public <V extends de.jtem.halfedge.Vertex<V, E, F>, E extends de.jtem.halfedge.Edge<V, E, F>, F extends de.jtem.halfedge.Face<V, E, F>> Double getE(
				E e, AdapterSet adapters) {

			if (!adapters.contains(Length.class, e.getClass(), Double.class)) {
				throw new RuntimeException(
						"Need adapter for length of edges to calculate cotan weights.");
			}

			double a = adapters.get(Length.class, e, Double.class);
			double b = adapters
					.get(Length.class, e.getNextEdge(), Double.class);
			double c = adapters.get(Length.class,
					e.getNextEdge().getNextEdge(), Double.class);

			if (!e.getNextEdge().getNextEdge().getNextEdge().equals(e)) {
				throw new RuntimeException("Face is not a triangle.");
			}

			double cosalpha = (b * b + c * c - a * a) / (2 * b * c);
			double cotanalpha = cosalpha / (Math.sqrt(1 - cosalpha * cosalpha));

			return cotanalpha;
		};

		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return Edge.class.isAssignableFrom(nodeClass);
		}

	}	

	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double[][] getHarmonicDifferentials(HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters) {
		// First make clear that we are working with a delaunay triangulation.
		DelaunayLengthAdapter la = Delaunay.constructDelaunay(hds, adapters);
		
		// Get the homology basis of the surface.
		V rootV = hds.getVertex(0);
		WeightAdapter<E> wa = new Search.DefaultWeightAdapter<E>();
		List<Set<E>> basis = HomologyUtility.getGeneratorPaths(rootV, wa);

		// For each cycle in the basis we get a corresponding harmonic
		// differential (closed but not exact, its integral differs by one on
		// the cycle).
		double[][] dh = new double[basis.size()][];

		// For each cycle construct the corresponding harmonic differential.
		adapters.add(la); // delaunay lengths
		for (int i = 0; i < basis.size(); i++) {
			dh[i] = getHarmonicDifferential(hds, adapters, basis.get(i));
		}
		adapters.remove(la);
		return dh;
	}
	
	
	private static enum EdgeStatus {
		EndsAtLeftSide,
		StartsAtLeftSide,
		EndsAtRightSide,
		StartsAtRightSide,
		liesOnTheCycle,
		NoConnection
	}
	
	
	
	/**
	 * Rotate e around v clockwise
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param e
	 * @param v
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> E getNextEdgeClockwise(E e, V v) {
		if (e.getStartVertex() == v) {
			return e.getPreviousEdge().getOppositeEdge();
		} if (e.getTargetVertex() == v) {
			return e.getOppositeEdge().getPreviousEdge();
		}
		throw new IllegalArgumentException("Edge does not contain vertex in getNextEdgeClockwise()");
	}
	
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> EdgeStatus getEdgeStatus(E e, Set<E> eCycle, Set<V> vCycle) {
		boolean outpointing = vCycle.contains(e.getStartVertex());
		if(eCycle.contains(e) || eCycle.contains(e.getOppositeEdge())) 
			return EdgeStatus.liesOnTheCycle;
		V v= outpointing ? e.getStartVertex(): e.getTargetVertex();
		E curr= getNextEdgeClockwise(e,v);
		while (curr == e){
			if(eCycle.contains(curr)) {
				if(outpointing) return EdgeStatus.StartsAtRightSide;
				else return EdgeStatus.EndsAtLeftSide;
			}
			if(eCycle.contains(curr.getOppositeEdge())) {
				if(outpointing) return EdgeStatus.StartsAtLeftSide;
				else return EdgeStatus.EndsAtRightSide;
			}
			curr=getNextEdgeClockwise(e,v);
		} 
		
		return EdgeStatus.NoConnection;
	}
	
	

	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double[] getHarmonicDifferential(
		HalfEdgeDataStructure<V, E, F> hds, 
		AdapterSet adapters,
		Set<E> cycle
	) {

		Map<Integer, Integer> tau = getIdentificationMap(hds, cycle);

		// Consider the surface hds being obtained by identifying two boundary
		// cycles of another surface M via the map above. On this surface we can
		// get a harmonic function with boundary conditions 0 and 1 constant on
		// the identified cycles.
		double[] h = getStandardHarmonicFunction(hds, adapters, cycle, tau);

		// The differential of h can be defined on hds.
		double[] dh = new double[hds.numEdges() / 2];
		Set<V> boundaryVertexSet = getVertexSet(hds, cycle);

		// build the differences for each edge
		for (E e : hds.getPositiveEdges()) {
			
			int k= adapters.get(EdgeIndex.class, e,Integer.class);
			EdgeStatus status= getEdgeStatus(e, cycle, boundaryVertexSet);
				
			switch (status) {
			case StartsAtRightSide:
				dh[k] = h[e.getTargetVertex().getIndex()]
							- h[tau.get(e.getStartVertex().getIndex())];
				break;
			case EndsAtRightSide:
				dh[k] = h[tau.get(e.getTargetVertex().getIndex())]
							- h[e.getStartVertex().getIndex()];
				break;
			default:
				dh[k]= h[e.getTargetVertex().getIndex()]
					- h[e.getStartVertex().getIndex()];
			break;
			}
		}

		return dh;
	}

	/**
	 * Returns a harmonic function, which is two valued on the boundary
	 * (constant 0 and 1 resp.).
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
	> double[] getStandardHarmonicFunction(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			Set<E> cycle, Map<Integer, Integer> tau) {
	
		double[] bc1 = new double[tau.size()];
		double[] bc2 = new double[tau.size()];
		for (int i = 0; i < bc2.length; i++) {
			bc2[i] = 1.;
		}
		return getHarmonicFunction(hds, adapters, cycle, tau, bc1, bc2);
		
	}

	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double[] getHarmonicFunction(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			Set<E> cycle, Map<Integer, Integer> tau,
			double[] boundaryCondition1, double[] boundaryCondition2) {

		int n = hds.numVertices() + cycle.size();
		CompRowMatrix laplaceop = (CompRowMatrix) getLaplacian(hds, adapters,
				cycle, tau);

		Set<V> vertexSet = getVertexSet(hds, cycle);

		Vector x = new DenseVector(n);
		Vector bc = new DenseVector(n);

		for (V v : vertexSet) {
			int i= v.getIndex();
			bc.set(i, boundaryCondition1[i]);
			bc.set(tau.get(i), boundaryCondition2[i]);
		}

		IterativeSolver solver = new GMRES(x);
		solver.getIterationMonitor().setIterationReporter(
				new OutputIterationReporter());

		try {
			solver.solve(laplaceop, bc, x);
		} catch (IterativeSolverNotConvergedException e) {
			System.err.println("Iterative solver failed to converge");
		}

		double[] H = new double[n];

		for (int i = 0; i < hds.numVertices(); i++) {
			H[i] = x.get(i);
		}

		return H;
	}

	/**
	 * Returns a map which maps the indices of the vertices v0,v1,... contained
	 * in the cycle to the integers n,n+1,... , where n is the number of
	 * vertices of the surface.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param cycle
	 * @return
	 */
	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> Map<Integer, Integer> getIdentificationMap(
			HalfEdgeDataStructure<V, E, F> hds, Set<E> cycle) {
		Map<Integer, Integer> map = new HashMap<Integer, Integer>(100);

		// Get the vertices contained in the cycle.
		Set<V> vertexSet = getVertexSet(hds, cycle);
		int n = hds.numVertices();

		// Build up the map.
		int i=0;
		for (V v:vertexSet) {
			map.put(v.getIndex(), n + i);
			i++;
		}
		return map;
	}

	/**
	 * Returns the indices of the vertices contained in the cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param cycle
	 * @return
	 */
	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> Set<V> getVertexSet(
			HalfEdgeDataStructure<V, E, F> hds, Set<E> cycle) {
		Set<V> vertexSet = new HashSet<V>(cycle.size());
		for (E e:cycle) {
			vertexSet.add(e.getStartVertex());
		}
		return vertexSet;
	}

	/**
	 * Returns the matrix of the cotan laplacian corresponding to the surface
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
	> Matrix getLaplacian(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			Set<E> cycle, Map<Integer, Integer> tau) {

		CotanAdapter ca = new CotanAdapter();
		adapters.add(ca);

		Set<V> boundaryVertexSet = getVertexSet(hds, cycle);

		DoubleMatrixFactory dmf = new DoubleMatrixFactory(hds.numVertices()+cycle.size(),
				hds.numVertices()+cycle.size());
		int i, j;
		double weight;

		// TODO: Implement correctly
		for (E e : hds.getEdges()) {
			int k= adapters.get(EdgeIndex.class, e,Integer.class);
			EdgeStatus status= getEdgeStatus(e, cycle, boundaryVertexSet);
				
//			switch (status) {
//			case StartsAtRightSide:
//				dh[k] = h[e.getTargetVertex().getIndex()]
//							- h[tau.get(e.getStartVertex().getIndex())];
//				break;
//			case EndsAtRightSide:
//				dh[k] = h[tau.get(e.getTargetVertex().getIndex())]
//							- h[e.getStartVertex().getIndex()];
//				break;
//			case liesOnTheCycle:
//				dmf.set(e.getStartVertex().getIndex(), e.getStartVertex().getIndex(), 1.);
//				dmf.set(tau.get(e.getStartVertex().getIndex()), e.getStartVertex().getIndex(), 1.);
//			default:
//				i = e.getStartVertex().getIndex();
//				j = e.getTargetVertex().getIndex();
//			break;
//			}
//			weight = getCotanWeight(e, adapters);
//			dmf.set(i, j, weight);
//			dmf.add(i, i, -weight);
			if (boundaryVertexSet.contains(e.getStartVertex())) {
				i = e.getStartVertex().getIndex();
				j = tau.get(e.getStartVertex().getIndex());
				dmf.set(i, i, 1.);
				dmf.set(j, j, 1.);
			} else {
				i = e.getStartVertex().getIndex();
				if (boundaryVertexSet.contains(e.getTargetVertex())
						&& !cycle.contains(e.getNextEdge()))
					j = tau.get(e.getTargetVertex().getIndex());
				else
					j = e.getTargetVertex().getIndex();
				weight = getCotanWeight(e, adapters);
				dmf.set(i, j, weight);
				dmf.add(i, i, -weight);
			}
		}
		
		dmf.update();
		adapters.remove(ca);
		return dmf.getMatrix();

	}

	private static double getCotanWeight(Edge<?, ?, ?> e, AdapterSet adapters) {
		return adapters.get(Weight.class, e, Double.class)
				+ adapters.get(Weight.class, e.getOppositeEdge(), Double.class);
	}

	/**
	 * This method shall return whether a given vertex lies on the boundary of
	 * the surface it belongs to.
	 * 
	 * @param v
	 * @return isBoundaryVertex
	 */
	private static boolean isBoundaryVertex(Vertex<?, ?, ?> v) {
		Edge<?, ?, ?> currEdge = v.getIncomingEdge();
		do {
			if (HalfEdgeUtilsExtra.isBoundaryEdge(currEdge))
				return true;
			currEdge = currEdge.getNextEdge().getOppositeEdge();
		} while (currEdge.equals(v.getIncomingEdge()));
		return false;
	}

}
