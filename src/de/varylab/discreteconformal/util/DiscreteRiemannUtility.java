package de.varylab.discreteconformal.util;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.BiCG;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import no.uib.cipr.matrix.sparse.ILU;
import no.uib.cipr.matrix.sparse.IterativeSolver;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import no.uib.cipr.matrix.sparse.OutputIterationReporter;
import no.uib.cipr.matrix.sparse.Preconditioner;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.EdgeIndex;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.Weight;
import de.jtem.halfedgetools.algorithm.triangulation.Delaunay;
import de.jtem.halfedgetools.algorithm.triangulation.DelaunayLengthAdapter;
import de.varylab.discreteconformal.util.Search.WeightAdapter;
import de.varylab.matrix.sparse.factory.DoubleMatrixFactory;

public class DiscreteRiemannUtility {

	/**
	 * Adapter returns the cotan weight of a given edge.
	 */
	@Weight
	private static class CotanAdapter extends AbstractAdapter<Double> {

		public CotanAdapter() {
			super(Double.class, true, false);
		}

		public <
			V extends Vertex<V, E, F>,
			E extends Edge<V, E, F>,
			F extends Face<V, E, F>
		> Double getE(E e, AdapterSet adapters) {
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
		public <
			N extends Node<?, ?, ?>
		> boolean canAccept(Class<N> nodeClass) {
			return Edge.class.isAssignableFrom(nodeClass);
		}

	}

	/**
	 * Calculates 2*g harmonic differentials on hds. The weight Adapter is used
	 * to find a basis of the homology that are short with respect to the weight
	 * given by wa.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param wa
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double[][] getHarmonicForms(
		HalfEdgeDataStructure<V, E, F> hds, 
		AdapterSet adapters,
		WeightAdapter<E> wa
	) {
		// First make clear that we are working with a delaunay triangulation.
		DelaunayLengthAdapter la = Delaunay.constructDelaunay(hds, adapters);

		// Get the homology basis of the surface.
		V rootV = hds.getVertex(0);
		List<Set<E>> basis = HomologyUtility.getGeneratorPaths(rootV, wa);

		// use the private method
		return getHarmonicForms(hds, basis, adapters, la, wa);
	}
	
	/**
	 * Returns a basis of holomorphic differentials on the surface hds. 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param wa
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double[][] getHolomorphicForms(
		HalfEdgeDataStructure<V, E, F> hds, 
		AdapterSet adapters,
		WeightAdapter<E> wa){
		
		// First make clear that we are working with a delaunay triangulation.
		DelaunayLengthAdapter la = Delaunay.constructDelaunay(hds, adapters);

		// Get the homology basis of the surface.
		V rootV = hds.getVertex(0);
		List<Set<E>> basis = HomologyUtility.getGeneratorPaths(rootV, wa);

		// use the private method
		return getHolomorphicForms(hds, basis, adapters, la, wa);
	}

	/**
	 * Returns a basis of holomorphic differentials on the surface, which is
	 * normalized with respect to the given homology basis. The surface has to
	 * be Delaunay triangulated.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param homologyBasis
	 * @param adapters
	 * @param la
	 * @param wa
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double[][] getHolomorphicForms(
		HalfEdgeDataStructure<V, E, F> delaunay,
		List<Set<E>> homologyBasis,
		AdapterSet adapters,
		DelaunayLengthAdapter la,
		WeightAdapter<E> wa){
		
		// TODO: implement me!

		double[][] dh = getHarmonicForms(delaunay, homologyBasis, adapters, la,
				wa);
		double[][] dhStar = getStarOfForms(delaunay, adapters, dh);

		List<Set<E>> acycles = getACycles(delaunay, homologyBasis);
		List<Set<E>> dualACycles = getDualPaths(delaunay, acycles);
		CompRowMatrix A = getCoefficientMatrix(dh, dhStar, acycles,
				dualACycles, adapters);

		DenseVector x = new DenseVector(dh.length);

		IterativeSolver solver = new BiCG(x);

		Preconditioner pre = new ILU(A);
		pre.setMatrix(A);
		solver.setPreconditioner(pre);

		solver.getIterationMonitor().setIterationReporter(
				new OutputIterationReporter());

		for (int i = 0; i < acycles.size(); i++) {
			DenseVector bc = new DenseVector(2 * acycles.size());
			bc.set(i, 1);
			try {
				solver.solve(A, bc, x);
			} catch (IterativeSolverNotConvergedException e) {
				System.err
						.println("Iterative solver failed to converge: Couldn't get harmonic function.");
			}

			// TODO: build linear combinations of the harmonic differentials
			// using the obtained coefficients, etc

		}

		return null;
	}
	
	/**
	 * Returns a paths in the dual surface which are homotopic to the given
	 * ones.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param cycles
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<Set<E>> getDualPaths(HalfEdgeDataStructure<V,E,F> hds, List<Set<E>> cycles){
		List<Set<E>> dualCycles= new java.util.Vector<Set<E>>();
		for (int i = 0; i < cycles.size(); i++) {
			dualCycles.add(getDualPath(hds,cycles.get(i)));
		}
		return dualCycles;
	}
	
	/**
	 * Returns a path in the dual surface which is homotopic to the given one.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param cycles
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<E> getDualPath(HalfEdgeDataStructure<V,E,F> hds, Set<E> cycles){
		// TODO: Implement me!
		return null;
	}
	
	/**
	 * Returns a matrix needed to build the normalized holomorphic differentials.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param dh
	 * @param dhStar
	 * @param acycles
	 * @param dualACycles
	 * @param adapters
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> CompRowMatrix getCoefficientMatrix(
			double[][] dh, double[][] dhStar, 
			List<Set<E>> acycles, List<Set<E>> dualACycles, 
			AdapterSet adapters){
		// TODO: Implement me!
		return null;
	}
	
	/**
	 * Returns g cycles of the homology basis representing the a cycles.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param homologyBasis
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> List<Set<E>> getACycles(HalfEdgeDataStructure<V,E,F> hds, List<Set<E>> homologyBasis){
		// TODO: Choose the a cycles and return them.
		return null;
	}
	
	/**
	 * Applies the Hodge star operator to an array of forms on the surface and
	 * returns the result.
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
	> double[][] getStarOfForms(
		HalfEdgeDataStructure<V, E, F> delaunay,
		AdapterSet adapters,
		double[][] forms){
		double[][] formsStar = new double[forms.length][];
		for (int i = 0; i < formsStar.length; i++) {
			formsStar[i] = getStarOfForm(delaunay, adapters, forms[i]);
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
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double[] getStarOfForm(
		HalfEdgeDataStructure<V, E, F> delaunay,
		AdapterSet adapters,
		double[] form){
		
		double[] formStar = new double[form.length];

		CotanAdapter ca = new CotanAdapter();

		double ratio;
		int id;

		for (E e : delaunay.getPositiveEdges()) {
			// TODO: Check! Is the cotan weight the ratio of the edge length and
			// its dual edge length?
			ratio = ca.getE(e, adapters);
			id = adapters.get(EdgeIndex.class, e, Integer.class);
			formStar[id] = ratio * form[id];
		}

		return formStar;
	}
	

	/**
	 * Calculates 2*g harmonic differentials on on a Delaunay triangulated
	 * surface hds corresponding to a specified homology basis. The weight
	 * Adapter is used to find a basis of the homology that are short with
	 * respect to the weight given by wa.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param homologyBasis
	 * @param adapters
	 * @param la
	 * @param wa
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double[][] getHarmonicForms(
		HalfEdgeDataStructure<V, E, F> delaunay, 
		List<Set<E>> homologyBasis,
		AdapterSet adapters, DelaunayLengthAdapter la,
		WeightAdapter<E> wa
	) {
		
		// For each cycle in the basis we get a corresponding harmonic
		// differential (closed but not exact, its integral differs by one on
		// the cycle).
		double[][] dh = new double[homologyBasis.size()][];

		// For each cycle construct the corresponding harmonic differential.
		adapters.add(la); // delaunay lengths
		for (int i = 0; i < homologyBasis.size(); i++) {
			dh[i] = getStandardHarmonicForm(delaunay, adapters,
					homologyBasis.get(i));
		}
		adapters.remove(la);
		return dh;
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
		}
		if (e.getTargetVertex() == v) {
			return e.getOppositeEdge().getPreviousEdge();
		}
		throw new IllegalArgumentException(
				"Edge does not contain vertex in getNextEdgeClockwise()");
	}
	
	/**
	 * There are several possibilities how an edge can be connected to the given
	 * cycle. If the surface is cut along the cycle, we can regard the cycle as
	 * two cycles actually - a left and a right one.
	 */
	private static enum EdgeStatus {
		endsAtLeftCycle,
		startsAtLeftCycle,
		endsAtRightCycle,
		startsAtRightCycle,
		liesOnLeftCycle,
		liesOnRightCycle,
		noConnection
	}

	/**
	 * Returns the edge's status, see EdgeStatus description.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param e
	 * @param edgeCycle
	 * @param vertexCycle
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> EdgeStatus getEdgeStatus(E e, Set<E> edgeCycle, Set<V> vertexCycle) {
		boolean outpointing = vertexCycle.contains(e.getStartVertex());
		if (edgeCycle.contains(e))
			return EdgeStatus.liesOnLeftCycle;
		if (edgeCycle.contains(e.getOppositeEdge()))
			return EdgeStatus.liesOnRightCycle;
		V v = outpointing ? e.getStartVertex() : e.getTargetVertex();
		E curr = getNextEdgeClockwise(e, v);
		while (curr != e) {
			if (edgeCycle.contains(curr)) {
				if (outpointing)
					return EdgeStatus.startsAtRightCycle;
				else
					return EdgeStatus.endsAtLeftCycle;
			}
			if (edgeCycle.contains(curr.getOppositeEdge())) {
				if (outpointing)
					return EdgeStatus.startsAtLeftCycle;
				else
					return EdgeStatus.endsAtRightCycle;
			}
			curr = getNextEdgeClockwise(curr, v);
		}

		return EdgeStatus.noConnection;
	}
	
	/**
	 * Returns the standard harmonic differential for the given cycle, i.e. the
	 * harmonic differential, whose integral along the cycle is 1. 
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param cycle
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double[] getStandardHarmonicForm(HalfEdgeDataStructure<V, E, F> hds, 
		AdapterSet adapters, Set<E> cycle) {

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
			case startsAtRightCycle:
				dh[k] = h[e.getTargetVertex().getIndex()]
							- h[tau.get(e.getStartVertex().getIndex())];
				break;
			case endsAtRightCycle:
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
	 * Returns a harmonic function, which is two valued on the cycle
	 * (constant 0 on the left and 1 on the right side).
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

	/**
	 * Returns a harmonic function on the surface cut along the given cycle. The
	 * cycle splits in a left and a right side, where boundary conditions can be
	 * specified. Tau is the identification map to obtain the surface given by
	 * the hds.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param cycle
	 * @param tau
	 * @param boundaryCondition1
	 *            (left side)
	 * @param boundaryCondition2
	 *            (right side)
	 * @return
	 */
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

		int bcIndex = 0;
		for (V v : vertexSet) {
			int i = v.getIndex();
			bc.set(i, boundaryCondition1[bcIndex]);
			bc.set(tau.get(i), boundaryCondition2[bcIndex]);
			bcIndex++;
		}

		IterativeSolver solver = new BiCG(x);
		
		Preconditioner pre= new ILU(laplaceop);
		pre.setMatrix(laplaceop);
		solver.setPreconditioner(pre);
		
		solver.getIterationMonitor().setIterationReporter(
				new OutputIterationReporter());

		try {
			solver.solve(laplaceop, bc, x);
		} catch (IterativeSolverNotConvergedException e) {
			System.err
					.println("Iterative solver failed to converge: Couldn't get harmonic function.");
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
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> Map<Integer, Integer> getIdentificationMap(
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
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<V> getVertexSet(
			HalfEdgeDataStructure<V, E, F> hds, Set<E> cycle) {
		Set<V> vertexSet = new HashSet<V>(cycle.size());
		for (E e:cycle) {
			vertexSet.add(e.getStartVertex());
		}
		return vertexSet;
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
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> Matrix getLaplacian(HalfEdgeDataStructure<V, E, F> hds, 
			AdapterSet adapters, Set<E> cycle, Map<Integer, Integer> tau) {

		Set<V> boundaryVertexSet = getVertexSet(hds, cycle);

		DoubleMatrixFactory dmf = new DoubleMatrixFactory(hds.numVertices()+cycle.size(),
				hds.numVertices()+cycle.size());
		
		EdgeStatus status;
		double weight;
		int i, j;

		CotanAdapter ca = new CotanAdapter();
		adapters.add(ca);
		
		// TODO: Check implementation!
		for (E e : hds.getEdges()) {
			status= getEdgeStatus(e, cycle, boundaryVertexSet);
			weight = getCotanWeight(e, adapters);
			switch (status) {
			case liesOnLeftCycle:
				i = e.getStartVertex().getIndex();
				dmf.set(i, i, 1.);
				break;
			case liesOnRightCycle:
				i = tau.get(e.getStartVertex().getIndex());
				dmf.set(i, i, 1.);
				break;
			case startsAtRightCycle:
				i = tau.get(e.getStartVertex().getIndex());
				j = e.getTargetVertex().getIndex();
				dmf.set(i, j, weight);
				dmf.add(i, i, -weight);
				break;
			case endsAtRightCycle:
				i = e.getStartVertex().getIndex();
				j = tau.get(e.getTargetVertex().getIndex());
				dmf.set(i, j, weight);
				dmf.add(i, i, -weight);
				break;
			case endsAtLeftCycle:
			default:
				i = e.getStartVertex().getIndex();
				j = e.getTargetVertex().getIndex();
				dmf.set(i, j, weight);
				dmf.add(i, i, -weight);
				break;
			}
		}
		
		adapters.remove(ca);
		
		dmf.update();
		return dmf.getMatrix();

	}

	private static double getCotanWeight(Edge<?, ?, ?> e, AdapterSet adapters) {
		return adapters.get(Weight.class, e, Double.class)
				+ adapters.get(Weight.class, e.getOppositeEdge(), Double.class);
	}

}
