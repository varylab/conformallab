package de.varylab.discreteconformal.util;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.solver.DefaultDoubleIterationMonitor;
import cern.colt.matrix.tdouble.algo.solver.DoubleBiCGstab;
import cern.colt.matrix.tdouble.algo.solver.DoubleIterationReporter;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

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
import de.jtem.halfedgetools.util.HalfEdgeUtilsExtra;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.util.Search.WeightAdapter;

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

	private static DoubleIterationReporter reporter = new DoubleIterationReporter() {
		
		@Override
		public void monitor(double arg0, DoubleMatrix1D arg1, int arg2) {
			monitor(arg0, arg2);
		}
		
		@Override
		public void monitor(double arg0, int arg1) {
			System.err.println("iteration = "+arg1+", value = "+arg0);
		}
	};
	
	private static double eps= 1E-10;
	
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

		// List<Set<E>> acycles = getACycles(hds, basis);
		//
		// for (Set<E> c: acycles) {
		// System.err.println(c);
		// }
		
		// use the private method
		return getHarmonicForms(hds, basis, adapters, la, wa);
	}

	/**
	 * Returns a basis of holomorphic differentials on the surface hds. The rows
	 * represent holomorphic forms. In the i-th column the value on the positive
	 * oriented edge with edgeIndex i is saved. The real part correspond to the
	 * edge the imaginary to its dual.
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
	> Complex[][] getHolomorphicForms(
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
	 * be Delaunay triangulated. The forms are saved as rows. The real part of
	 * the value in the i-th row corresponds to the value of the differential on
	 * the positive oriented edge with edgeIndex i, the imaginary part to the
	 * value of its dual edge.
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
	> Complex[][] getHolomorphicForms(
		HalfEdgeDataStructure<V, E, F> delaunay,
		List<Set<E>> homologyBasis,
		AdapterSet adapters,
		DelaunayLengthAdapter la,
		WeightAdapter<E> wa){

		// the forms are defined on the positive oriented edges
		int numEdges = delaunay.numEdges()/2;

		// Get the harmonic differentials on the surface and its dual. The
		// format of the matrices is 2g*numEdges
		double[][] dh = getHarmonicForms(delaunay, homologyBasis, adapters, la,
				wa);
		double[][] dhStar = getDualForms(delaunay, adapters, dh);

		// needed to multiply colt matrices
		DenseDoubleAlgebra dalgebra = new DenseDoubleAlgebra();

		// We want to have matrices each column of which corresponds to one
		// harmonic differential. So we can build linear combinations of the
		// harmonic differentials simply by matrix multiplication. TODO: Check,
		// whether, transpose is perhaps too expensive?
		DoubleMatrix2D dhTransposed = dalgebra.transpose(DoubleFactory2D.dense
				.make(dh));
		DoubleMatrix2D dhStarTransposed = dalgebra
				.transpose(DoubleFactory2D.dense.make(dhStar));

		// to normalize the differentials we need the a periods and its dual
		// cycles
		List<Set<E>> acycles = getACycles(delaunay, homologyBasis);
		List<Set<E>> dualACycles = getDualPaths(delaunay, acycles);

		// g is simply the genus of the surface
		int g = acycles.size();

		// an array to save the forms: The real part corresponds to the
		// differential on the primal mesh and the imaginary part to the
		// differentials on the dual mesh. 
		Complex[][] OMEGA= new Complex[g][numEdges];
		
		// The holomorphic forms omega are linear combinations of 2g harmonic
		// forms dh_j plus i times its duals (star)dh_j: omega= sum(x_j*(dh_j+
		// i*(star)dh_j)). They are normalized iff they satisfy a linear system.
		// A is the coefficient matrix.
		DoubleMatrix2D A = getCoefficientMatrix(dh, dhStar, acycles,
				dualACycles, adapters);

		// The number of harmonic forms on the surface is 2g, so we need 2g
		// coefficients.
		DoubleMatrix1D x = DoubleFactory1D.dense.make(2 * g);

		// Iterative solver
		DoubleBiCGstab solver = new DoubleBiCGstab(x);
		DefaultDoubleIterationMonitor monitor= new DefaultDoubleIterationMonitor();
		
		// configure monitor
		monitor.setMaxIterations(100000);
		monitor.setAbsoluteTolerance(eps);
		monitor.setRelativeTolerance(eps);
		monitor.setIterationReporter(reporter);
		
		solver.setIterationMonitor(monitor);

		// For each a-cycle find the holomorphic form which is 1 along this
		// cycle, i.e. gives one for the primal cycle and 0 for its dual cycle. 
		for (int i = 0; i < g; i++) {
			
			// set up the conditions 
			DoubleMatrix1D bc = DoubleFactory1D.dense.make(2 * g);
			bc.set(i, 1);
			
			// solve the system
			try {
				solver.solve(A, bc, x);
			} catch (IterativeSolverDoubleNotConvergedException e) {
				System.err
						.println("Iterative solver failed to converge: Couldn't get holomorphic form.");
				e.printStackTrace();
			}

			// Build linear combinations of the harmonic differentials
			// using the obtained coefficients. 
			DoubleMatrix1D omega= dalgebra.mult(dhTransposed, x);
			DoubleMatrix1D omegaStar= dalgebra.mult(dhStarTransposed, x);
			
			for (int j = 0; j < numEdges; j++) {
				OMEGA[i][j]= new Complex(omega.get(j),omegaStar.get(j));
			}
		}

		return OMEGA;
	}
	
	/**
	 * Returns paths in the dual surface which are homotopic to the given
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
		// for each cycle
		for (int i = 0; i < cycles.size(); i++) {
			dualCycles.add(getDualPath(hds, cycles.get(i)));
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
	 * @param cycle
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<E> getDualPath(HalfEdgeDataStructure<V,E,F> hds, Set<E> cycle){
		Set<E> dualPath = new HashSet<E>();
		Set<V> vertices = getVertexSet(cycle);
		for (V v : vertices) {
			List<E> star = HalfEdgeUtilsExtra.getEdgeStar(v);
			for (E e : star) {
				EdgeStatus status = getEdgeStatus(e, cycle, vertices);
				if (status == EdgeStatus.endsAtLeftCycle
						|| status == EdgeStatus.startsAtLeftCycle)
					dualPath.add(e);
			}
		}
		return dualPath;
	}

	/**
	 * Returns a matrix A needed to build the normalized holomorphic
	 * differentials. The matrix has format 2g*2g. There are 2g harmonic
	 * differentials dH_j= dh_j+i*(star)dh_j. The first g rows are filled with
	 * its a-periods, the second g rows with its values on the a-cycles of dual
	 * surface. That means in the upper g rows of Ax stands the real part and
	 * in the lower the imaginary part of the a-periods of the holomorphic
	 * differential omega=sum(x_j*dH_j).
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param dh
	 * @param dhStar
	 * @param aCycles
	 * @param dualACycles
	 * @param adapters
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix2D getCoefficientMatrix(
			double[][] dh, double[][] dhStar, 
			List<Set<E>> aCycles, List<Set<E>> dualACycles, 
			AdapterSet adapters){
		
		// The genus equals the number of a-cycles
		int g= aCycles.size();

		// The matrix has to have the format 2g*2g. TODO: Check whether the
		// matrix should better be defined dense. 
		DoubleMatrix2D M= DoubleFactory2D.sparse.make(2*g, 2*g);
		
		// for each cycle and each differential
		for (int i = 0; i < g; i++) {
			for (int j = 0; j < 2 * g; j++) {
				// fill the upper g rows with a-periods of the primal mesh
				M.set(i, j, integrateFormOverCycle(dh[j], aCycles.get(i),
						adapters));
				// and the lower g rows with the a-periods of the dual mesh
				M.set(i + g, j, integrateFormOverCycle(dhStar[j], dualACycles
						.get(i), adapters));
			}
		}
		
		return M;
	}

	/**
	 * Returns the integral of the differential over the given cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param form
	 * @param cycle
	 * @param adapters
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double integrateFormOverCycle(
			double[] form, 
			Set<E> cycle,
			AdapterSet adapters){
		// init with zero
		double integral = 0;
		int id;
		// for each edge of the cycle
		for (E e : cycle) {
			// get global index
			id = adapters.get(EdgeIndex.class, e, Integer.class);
			if (e.isPositive()) // if the edge is positive add
				integral += form[id];
			else
				// if negative subtract the value of the form on the edge
				integral -= form[id];
		}
		return integral;
	}

	/**
	 * Returns g cycles of the homology basis representing the a cycles. TODO:
	 * Needs a canonical basis? It's not clear whether this really works in
	 * general. Why there shall be always g such cycles for an arbitrary basis.
	 * Would be good to visualize the returned cycles on the surface.
	 * 
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
		
		// the genus has to be one at least
		if (homologyBasis.size() < 1)
			return null;

		// a basis contains 2g cycles
		int g = homologyBasis.size() / 2;

		// g non-intersecting cycles shall be labeled as a-cycles
		List<Set<E>> aCycles = new Vector<Set<E>>(g);
		aCycles.add(homologyBasis.get(0));

		// try to find g cycles, which does not intersect
		for (int i = 1; i < g; i++) {
			// try to find a new cycle with intersection number 0
			for (Set<E> c : homologyBasis) {
				// if the cycle is already in the list, ignore it
				if (aCycles.contains(c))
					continue;
				// say c is non-intersecting
				boolean nonintersecting = true;
				// check the intersection number with each cycle already
				// collected
				for (Set<E> a : aCycles) {
					// if the intersection number is not equal to zero, c is
					// intersecting and has to be left out
					if (getIntersectionNumber(c, a) != 0) {
						nonintersecting = false;
						break;
					}
				}
				// if c is not intersecting any of the already collected cycles
				// add it to the list
				if (nonintersecting) {
					aCycles.add(c);
					break;
				}
			}
		}

		if (aCycles.size() != g)
			throw new RuntimeException(
					"Couldn't find g non-intersecting cycles.");

		return aCycles;
	}
	
	/**
	 * Returns the intersection number of two oriented cycles.
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param cycle1
	 * @param cycle2
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> int getIntersectionNumber(Set<E> cycle1, Set<E> cycle2) {
		// get the number of positive intersections of the first cycle with the
		// second
		int numPositiveIntersections = countEdgesWithStatus(cycle1, cycle2,
				EdgeStatus.startsAtLeftCycle);
		// get the number of negative intersections of the first cycle with the
		// second
		int numNegativeIntersections = countEdgesWithStatus(cycle1, cycle2,
				EdgeStatus.endsAtLeftCycle);
		return numPositiveIntersections - numNegativeIntersections;
	}

	/**
	 * Count the edges in a given set having a specified status with respect to
	 * a given cycle.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param cycle
	 * @param set
	 * @param status
	 * @return
	 */
	private static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> int countEdgesWithStatus(Set<E> cycle, Set<E> set, EdgeStatus status){
		int count = 0;
		Set<V> vertexSet= getVertexSet(cycle);
		for(E e:set){
			if(getEdgeStatus(e, cycle,vertexSet)== status)
				count++;
		}
		return count;
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
	> double[][] getDualForms(
		HalfEdgeDataStructure<V, E, F> delaunay,
		AdapterSet adapters,
		double[][] forms){
		double[][] formsStar = new double[forms.length][];
		for (int i = 0; i < formsStar.length; i++) {
			formsStar[i] = getDualForm(delaunay, adapters, forms[i]);
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
	> double[] getDualForm(
		HalfEdgeDataStructure<V, E, F> delaunay,
		AdapterSet adapters,
		double[] form){
		
		double[] formStar = new double[form.length];

		double ratio;
		int id;

		for (E e : delaunay.getPositiveEdges()) {
			// TODO: Check! Is the cotan weight the ratio of the edge length and
			// its dual edge length?
			ratio = getCotanWeight(e, adapters);
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
		Set<V> boundaryVertexSet = getVertexSet(cycle);

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
				dh[k] = h[e.getTargetVertex().getIndex()]
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
			bc2[i] = 10.;
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
		DoubleMatrix2D laplaceop = getLaplacian(hds, adapters, cycle, tau);

		Set<V> vertexSet = getVertexSet(cycle);

		DoubleMatrix1D x = DoubleFactory1D.dense.make(n);

		double[] bcond = new double[n];
		int bcIndex = 0;
		for (V v : vertexSet) {
			int id = v.getIndex();
			bcond[id] = boundaryCondition1[bcIndex];
			bcond[tau.get(id)] = boundaryCondition2[bcIndex];
			bcIndex++;
		}

		DoubleMatrix1D b = DoubleFactory1D.dense.make(bcond);

		DoubleBiCGstab solver = new DoubleBiCGstab(x);
		DefaultDoubleIterationMonitor monitor= new DefaultDoubleIterationMonitor();
		
		// configure monitor
		monitor.setMaxIterations(100000);
		monitor.setAbsoluteTolerance(eps);
		monitor.setRelativeTolerance(eps);
		monitor.setIterationReporter(reporter);
		
		solver.setIterationMonitor(monitor);

		try {
			solver.solve(laplaceop, b, x);
		} catch (IterativeSolverDoubleNotConvergedException e) {
			System.err
					.println("Iterative solver failed to converge: Couldn't get harmonic function.");
			e.printStackTrace();
		}
		System.err.println();

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
		Set<V> vertexSet = getVertexSet(cycle);
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
	> Set<V> getVertexSet(Set<E> cycle) {
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
	> DoubleMatrix2D getLaplacian(HalfEdgeDataStructure<V, E, F> hds, 
			AdapterSet adapters, Set<E> cycle, Map<Integer, Integer> tau) {

		Set<V> boundaryVertexSet = getVertexSet(cycle);

		DoubleMatrix2D M= DoubleFactory2D.sparse.make(hds.numVertices()+cycle.size(),
				hds.numVertices()+cycle.size());
		
		EdgeStatus status;
		double weight;
		int i, j;

		CotanAdapter ca = new CotanAdapter();
		adapters.add(ca);
		
		for (E e : hds.getEdges()) {
			status= getEdgeStatus(e, cycle, boundaryVertexSet);
			weight = getCotanWeight(e, adapters);
			switch (status) {
			case liesOnLeftCycle:
				i = e.getStartVertex().getIndex();
				M.set(i, i, 1.);
				break;
			case liesOnRightCycle:
				i = tau.get(e.getStartVertex().getIndex());
				M.set(i, i, 1.);
				break;
			case startsAtRightCycle:
				i = tau.get(e.getStartVertex().getIndex());
				j = e.getTargetVertex().getIndex();
				M.set(i, j, weight);
				M.set(i, i, M.get(i, i) - weight);
				break;
			case endsAtRightCycle:
				i = e.getStartVertex().getIndex();
				j = tau.get(e.getTargetVertex().getIndex());
				M.set(i, j, weight);
				M.set(i, i, M.get(i, i) - weight);
				break;
			case endsAtLeftCycle:
			default:
				i = e.getStartVertex().getIndex();
				j = e.getTargetVertex().getIndex();
				M.set(i, j, weight);
				M.set(i, i, M.get(i, i) - weight);
				break;
			}
		}
		
		adapters.remove(ca);
		return M;

	}

	private static double getCotanWeight(Edge<?, ?, ?> e, AdapterSet adapters) {
		return .5 * (adapters.get(Weight.class, e, Double.class) + adapters
				.get(Weight.class, e.getOppositeEdge(), Double.class));
	}

}
