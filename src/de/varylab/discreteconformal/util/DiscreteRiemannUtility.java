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
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.Weight;
import de.jtem.halfedgetools.adapter.type.generic.Position3d;
import de.jtem.halfedgetools.util.HalfEdgeUtilsExtra;
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

	/**
	 * Returns the R3 valued harmonic function fixing the boundary of the
	 * current immersion. (Only to test the cotan laplace.)
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param heds
	 * @param adapters
	 * @return
	 */
	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> double[][] getHarmonic(
			HalfEdgeDataStructure<V, E, F> heds, AdapterSet adapters) {
		int n = heds.numVertices();
		CompRowMatrix laplaceop = (CompRowMatrix) getLaplacian(heds, adapters);
		// PrintUtils.print(laplaceop);
		// PrintUtils.print(getLaplacian(heds, adapters));
		BlockMatrixFactory<Double> factory = new BlockMatrixFactory<Double>(n,
				n, 3, 3) {
			@Override
			protected double[][] getRealMatrixRepresentation(Double entry) {
				return new double[][] { { entry, 0, 0 }, { 0, entry, 0 },
						{ 0, 0, entry } };
			}

			@Override
			protected Double add(Double current, Double toAdd) {
				return current + toAdd;
			}
		};
		for (MatrixEntry e : laplaceop) {
			// System.err.println(e.row());
			factory.set(e.row(), e.column(), e.get());
		}

		factory.update();
		// PrintUtils.print(factory.getMatrix());
		// for (MatrixEntry e : factory.getMatrix()) {
		// System.err.println(e.row());
		// factory.set(e.row(), e.column(), e.get());
		// }
		CompRowMatrix A = (CompRowMatrix) factory.getMatrix();
		Vector X = new DenseVector(3 * n);
		Vector B = new DenseVector(3 * n);

		for (int i = 0; i < heds.numVertices(); i++) {
			if (isBoundaryVertex(heds.getVertex(i))) {
				double[] p = adapters.get(Position3d.class, heds.getVertex(i),
						double[].class);
				for (int k = 0; k < 3; k++)
					B.set(3 * i + k, p[k]);
			}
		}

		IterativeSolver solver = new GMRES(X);
		solver.getIterationMonitor().setIterationReporter(
				new OutputIterationReporter());

		try {
			solver.solve(A, B, X);
		} catch (IterativeSolverNotConvergedException e) {
			System.err.println("Iterative solver failed to converge");
		}

		double[][] F = new double[n][3];

		for (int i = 0; i < heds.numVertices(); i++) {
			for (int k = 0; k < 3; k++)
				F[i][k] = X.get(3 * i + k);
		}

		return F;
	}

	/**
	 * Returns the Delaunay triangulation of the given surface. That means all
	 * the cotan weights will be positive after that.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @return
	 */
	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> HalfEdgeDataStructure<V, E, F> getDelaunayTriangulation(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters) {
		// TODO: Implement it!
		return null;
	}

	/**
	 * Returns a basis of homology. It is represented by integer vectors that
	 * refer to the positive edges of the surface.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @return
	 */
	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> List<Set<E>> getHomologyBasis(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters) {
		// TODO: Implement it!
		return null;
	}

	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> double[][] getHarmonicDifferentials(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters) {

		// First make clear that we are working with a delaunay triangulation.
		HalfEdgeDataStructure<V, E, F> delaunay = getDelaunayTriangulation(hds,
				adapters);

		// Get the homology basis of the surface.
		List<Set<E>> basis = getHomologyBasis(delaunay, adapters);

		// For each cycle in the basis we get a corresponding harmonic
		// differential (closed but not exact, its integral differs by one on
		// the cycle).
		double[][] dh = new double[basis.size()][];

		// For each cycle construct the corresponding harmonic differential.
		for (int i = 0; i < basis.size(); i++) {
			dh[i] = getHarmonicDifferential(delaunay, adapters, basis.get(i));
		}

		return dh;
	}

	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> double[] getHarmonicDifferential(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			Set<E> cycle) {

		Map<Integer, Integer> tau = getIdentificationMap(hds, cycle);

		// Consider the surface hds beeing obtained by identifiying two boundary
		// cycles of another surface M via the map above. On this surface we can
		// get a harmonic function with boundary conditions 0 and 1 constant on
		// the identified cycles.
		double[] h = getStandardHarmonicFunction(hds, adapters, cycle, tau);

		// The differential of h can be defined on hds.
		double[] dh = new double[hds.numEdges() / 2];
		Set<V> boundaryVertexSet = getVertexSet(hds, cycle);

		// To create an int variable k, which run with, is not a solution. TODO:
		// Get out, how to write down a one form.
		int k = 0;

		// TODO: get the differential from the harmonic function.
		// ...
		for (E e : hds.getPositiveEdges()) {
			if (cycle.contains(e))
				// Function is constant on the cycle, i.e. dh(e) = 0.
				continue;
			else {
				// If the edge is not in the cycle, we still have to care about
				// edges, which touch the cycle, i.e. its start or target vertex
				// belong to the cycle.
				// On one side the values on the cycle differ.
				if (boundaryVertexSet.contains(e.getStartVertex())) {
					if (cycle.contains(e.getPreviousEdge()))
						dh[k] = h[e.getTargetVertex().getIndex()]
								- h[e.getStartVertex().getIndex()];
					else
						dh[k] = h[e.getTargetVertex().getIndex()]
								- h[tau.get(e.getStartVertex().getIndex())];
				} else if (boundaryVertexSet.contains(e.getTargetVertex())) {
					if (cycle.contains(e.getNextEdge()))
						dh[k] = h[e.getTargetVertex().getIndex()]
								- h[e.getStartVertex().getIndex()];
					else
						dh[k] = h[tau.get(e.getTargetVertex().getIndex())]
								- h[e.getStartVertex().getIndex()];
				} else {
					dh[k] = h[e.getTargetVertex().getIndex()]
							- h[e.getStartVertex().getIndex()];
				}
			}
			k++;
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
	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> double[] getStandardHarmonicFunction(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			Set<E> cycle, Map<Integer, Integer> tau) {
		double[] bc1 = new double[tau.size()];
		double[] bc2 = new double[tau.size()];
		for (int i = 0; i < bc2.length; i++) {
			bc2[i] = 1.;
		}
		return getHarmonicFunction(hds, adapters, cycle, tau, bc1, bc2);
	}

	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> double[] getHarmonicFunction(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			Set<E> cycle, Map<Integer, Integer> tau,
			double[] boundaryCondition1, double[] boundaryCondition2) {

		int n = hds.numVertices() + cycle.size();
		CompRowMatrix laplaceop = (CompRowMatrix) getLaplacian(hds, adapters,
				cycle, tau);

		Set<V> vertexSet = getVertexSet(hds, cycle);

		Vector x = new DenseVector(n);
		Vector bc = new DenseVector(n);

		for (int i = 0; i < hds.numVertices(); i++) {
			if (vertexSet.contains(hds.getVertex(i))) {
				bc.set(i, boundaryCondition1[i]);
				bc.set(tau.get(i), boundaryCondition2[i]);
			}
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
	 * Returns the Matrix of the standard cotan laplacian.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param heds
	 * @param adapters
	 * @return
	 */
	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> Matrix getLaplacian(
			HalfEdgeDataStructure<V, E, F> heds, AdapterSet adapters) {

		CotanAdapter ca = new CotanAdapter();
		adapters.add(ca);

		DoubleMatrixFactory dmf = new DoubleMatrixFactory(heds.numVertices(),
				heds.numVertices());
		int i, j;
		double weight;

		for (E e : heds.getEdges()) {
			if (isBoundaryVertex(e.getStartVertex()))
				dmf.set(e.getStartVertex().getIndex(), e.getStartVertex()
						.getIndex(), 1.);
			else {
				i = e.getStartVertex().getIndex();
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
	public static <V extends Vertex<V, E, F>, E extends Edge<V, E, F>, F extends Face<V, E, F>> Matrix getLaplacian(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			Set<E> cycle, Map<Integer, Integer> tau) {

		CotanAdapter ca = new CotanAdapter();
		adapters.add(ca);

		Set<V> boundaryVertexSet = getVertexSet(hds, cycle);

		DoubleMatrixFactory dmf = new DoubleMatrixFactory(hds.numVertices(),
				hds.numVertices());
		int i, j;
		double weight;

		// TODO: Check this again!
		for (E e : hds.getEdges()) {
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
