package de.varylab.discreteconformal.util;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import cern.colt.matrix.Norm;
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.solver.DefaultDoubleIterationMonitor;
import cern.colt.matrix.tdouble.algo.solver.DoubleGMRES;
import cern.colt.matrix.tdouble.algo.solver.DoubleIterationReporter;
import cern.colt.matrix.tdouble.algo.solver.DoubleIterativeSolver;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.EdgeIndex;
import de.jtem.halfedgetools.algorithm.triangulation.Delaunay;
import de.jtem.halfedgetools.algorithm.triangulation.MappedLengthAdapter;
import de.jtem.halfedgetools.util.HalfEdgeUtilsExtra;
import de.varylab.discreteconformal.util.EdgeUtility.EdgeStatus;
import de.varylab.discreteconformal.util.Search.WeightAdapter;


/**
 * Class to calculate harmonic differentials corresponding to a cycle of the
 * homology basis.
 * 
 * By convention all the forms are saved as rows and all the cycles as columns.
 * (For private methods we always use colt matrices.)
 * 
 * @author knoeppel
 * 
 */
public class DiscreteHarmonicFormUtility {

	private static DoubleIterationReporter reporter = new DoubleIterationReporter() {
		
		@Override
		public void monitor(double arg0, DoubleMatrix1D arg1, int arg2) {
			monitor(arg0, arg2);
		}

		@Override
		public void monitor(double arg0, int arg1) {
			if (arg1 % 100 == 0)
				System.err.println("iteration = " + arg1 + ", value = " + arg0);
		}
	};
	
	private static DenseDoubleAlgebra dalgebra = new DenseDoubleAlgebra();
	private static double eps = 1E-10;
	private static int maxIterations= 100000000;
	
	/**
	 * Calculates 2*g harmonic differentials on hds. The weight Adapter is used
	 * to find a basis of the homology that are short with respect to the weight
	 * given by wa. The forms are stored in the rows.
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
	> double[][] getHarmonicFormsOnPrimalMesh(
		HalfEdgeDataStructure<V, E, F> hds, 
		AdapterSet adapters,
		WeightAdapter<E> wa
	) {
		// First make clear that we are working with a delaunay triangulation.
		MappedLengthAdapter la = Delaunay.constructDelaunay(hds, adapters);

		// Get the homology basis of the surface.
		V rootV = hds.getVertex(0);
		List<List<E>> basis = CanonicalBasisUtility.getCanonicalHomologyBasis(
				rootV, adapters, wa);

		DoubleMatrix2D dh= getHarmonicFormsOfPrimalMesh(hds, basis, adapters, la, wa);

		// print(dalgebra.mult(dh, cyclesToMatrix(adapters, hds,basis)), 4);
		
		// use the private method
		return dh.toArray();
	}

	/**
	 * Calculates 2*g harmonic differentials on the dual of hds. The weight
	 * Adapter is used to find a basis of the homology that are short with
	 * respect to the weight given by wa. The forms are stored in the rows.
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
	> double[][] getHarmonicFormsOnDualMesh(
		HalfEdgeDataStructure<V, E, F> hds, 
		AdapterSet adapters,
		WeightAdapter<E> wa
	) {
		// First make clear that we are working with a delaunay triangulation.
		MappedLengthAdapter la = Delaunay.constructDelaunay(hds, adapters);

		// Get the homology basis of the surface.
		V rootV = hds.getVertex(0);
		List<List<E>> dualbasis = DualityUtility.getDualPaths(hds,
				CanonicalBasisUtility.getCanonicalHomologyBasis(rootV,
						adapters, wa));

		DoubleMatrix2D dh = getHarmonicFormsOfDualMesh(hds, dualbasis, adapters, la, wa);

		// print(dalgebra.mult(dh, cyclesToMatrix(adapters, hds,basis)), 4);
		
		// use the private method
		return dh.toArray();
	}

	/**
	 * Returns a matrix A needed to build the normalized holomorphic
	 * differentials. The matrix has format 2g*2g. There are 2g harmonic
	 * differentials dH_j= dh1_j+i*dh2_j. The first g rows are filled with
	 * its a1-periods, the second g rows with its values on the a2-cycles. 
	 * That means in the upper g rows of Ax stands the real part and in
	 * the lower the imaginary part of the a-periods of the holomorphic
	 * differential omega=sum(x_j*dH_j).
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param DH1
	 * @param DH2
	 * @param A1
	 * @param A2
	 * @param adapters
	 * @return
	 */
	public static DoubleMatrix2D getHarmonicPeriodsMatrix(
			DoubleMatrix2D DH1, DoubleMatrix2D DH2, DoubleMatrix2D A1,
			DoubleMatrix2D A2) {

		// The genus equals the number of a-cycles equals the number of columns of a
		int g = A1.columns();

		// get a-period-matrix of the harmonic differentials DH
		DoubleMatrix2D aPeriods1 = dalgebra.mult(DH1,A1);
		DoubleMatrix2D aPeriods2 = dalgebra.mult(DH2,A2);

		// The matrix has to have the format 2g*2g. 
		DoubleMatrix2D M = DoubleFactory2D.sparse.make(2 * g, 2 * g);

		// for each cycle and each differential
		for (int i = 0; i < g; i++) {
			for (int j = 0; j < 2 * g; j++) {
				// fill the upper g rows with a-periods of the primal mesh
				M.set(i, j, aPeriods1.get(j, i));
				// and the lower g rows with the a-periods of the dual mesh
				M.set(i + g, j, aPeriods2.get(j, i));
			}
		}

		// print(M, 4);

		return M;
	}
	
	/**
	 * Returns how harmonic a differential is. //TODO: find the right criterion,
	 * which is also easy to evaluate
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param hds
	 * @param adapters
	 * @param form
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double howHarmonic(HalfEdgeDataStructure<V,E,F> hds, AdapterSet adapters, DoubleMatrix1D form){
		double abssum= 0;
		int id;
		for (V v: hds.getVertices()) {
			List<E> star= HalfEdgeUtilsExtra.getEdgeStar(v);
			double val= 0;
			for (E e: star) {
				id= adapters.get(EdgeIndex.class, e,Integer.class);
				if(e.isPositive())
					val+= form.get(id);
				else
					val-= form.get(id);
			}
			abssum+=Math.abs(val);
		}
		return abssum;
	}

	/**
	 * Calculates 2*g harmonic differentials on on a Delaunay triangulated
	 * surface hds corresponding to a specified homology basis. The weight
	 * Adapter is used to find a basis of the homology that are short with
	 * respect to the weight given by wa. One differential corresponds to the
	 * entries in a row.
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
	> DoubleMatrix2D getHarmonicFormsOfPrimalMesh(
				HalfEdgeDataStructure<V, E, F> delaunay,
				List<List<E>> homologyBasis, AdapterSet adapters,
				MappedLengthAdapter la, WeightAdapter<E> wa) {
		
		// For each cycle in the basis we get a corresponding harmonic
		// differential (closed but not exact, its integral differs by one on
		// the cycle).
		DoubleMatrix2D dh = DoubleFactory2D.sparse.make(homologyBasis.size(),
				delaunay.numEdges() / 2);

		// For each cycle construct the corresponding harmonic differential.
		adapters.add(la); // delaunay lengths
		
		DoubleMatrix1D form;
		for (int i = 0; i < homologyBasis.size(); i++) {
			form= getStandardHarmonicFormOnPrimalMesh(delaunay, adapters,
					homologyBasis.get(i));
			for (int j = 0; j < form.size(); j++) {
				if (form.get(j) != 0)
					dh.set(i, j, form.get(j));
			}
		}
		adapters.remove(la);
		
		return dh;
	}

	/**
	 * Calculates 2*g harmonic differentials on on the dual a Delaunay
	 * triangulated surface hds corresponding to a specified homology basis. The
	 * weight Adapter is used to find a basis of the homology that are short
	 * with respect to the weight given by wa. One differential corresponds to
	 * the entries in a row.
	 * 
	 * @param <V>
	 * @param <E>
	 * @param <F>
	 * @param delaunay
	 * @param dualHomologyBasis
	 * @param adapters
	 * @param la
	 * @param wa
	 * @return
	 */
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> DoubleMatrix2D getHarmonicFormsOfDualMesh(
				HalfEdgeDataStructure<V, E, F> delaunay,
				List<List<E>> dualHomologyBasis, AdapterSet adapters,
				MappedLengthAdapter la, WeightAdapter<E> wa) {
		
		// For each cycle in the basis we get a corresponding harmonic
		// differential (closed but not exact, its integral differs by one on
		// the cycle).
		DoubleMatrix2D dh = DoubleFactory2D.sparse.make(dualHomologyBasis.size(),
				delaunay.numEdges() / 2);

		// For each cycle construct the corresponding harmonic differential.
		adapters.add(la); // delaunay lengths
		
		DoubleMatrix1D form;
		for (int i = 0; i < dualHomologyBasis.size(); i++) {
			form= getStandardHarmonicFormOnDualMesh(delaunay, adapters,
					dualHomologyBasis.get(i));
			for (int j = 0; j < form.size(); j++) {
				if (form.get(j) != 0)
					dh.set(i, j, form.get(j));
			}
		}
		
		DoubleMatrix2D Cycles= EdgeUtility.cyclesToMatrix(adapters, delaunay, dualHomologyBasis);
		SimpleMatrixPrintUtility.print(dalgebra.mult(dh, Cycles),4);
		
		adapters.remove(la);
		
		return dh;
	}
	
	/**
	 * Returns the standard harmonic differential on the primal mesh for the
	 * given cycle, i.e. the harmonic differential, whose integral along the
	 * cycle is 1.
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
	> DoubleMatrix1D getStandardHarmonicFormOnPrimalMesh(
				HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
				List<E> cycle) {

		Map<Integer, Integer> tau = getIdentificationMapOnPrimalMesh(hds, cycle);

		// Consider the surface hds being obtained by identifying two boundary
		// cycles of another surface M via the map above. On this surface we can
		// get a harmonic function with boundary conditions 0 and 1 constant on
		// the identified cycles.
		DoubleMatrix1D h = getStandardHarmonicFunctionOnPrimalMesh(hds, adapters, cycle,
				tau);

		// The differential of h can be defined on hds.
		DoubleMatrix1D dh = DoubleFactory1D.dense.make(hds.numEdges() / 2);
		Set<V> boundaryVertexSet = EdgeUtility.getPrimalVertexSet(cycle);

		// build the differences for each edge
		for (E e : hds.getPositiveEdges()) {

			int k = adapters.get(EdgeIndex.class, e, Integer.class);
			EdgeStatus status = EdgeUtility.getPrimalEdgeStatus(e, cycle,boundaryVertexSet);

			switch (status) {
			case startsAtRightCycle:
				dh.set(k,
						h.get(e.getTargetVertex().getIndex())
								- h.get(tau.get(e.getStartVertex().getIndex())));
				break;
			case endsAtRightCycle:
				dh.set(k,
						h.get(tau.get(e.getTargetVertex().getIndex()))
								- h.get(e.getStartVertex().getIndex()));
				break;
			default:
				dh.set(k,
						h.get(e.getTargetVertex().getIndex())
								- h.get(e.getStartVertex().getIndex()));
				break;
			}
		}

		return dh;
	}

	/**
	 * Returns the standard harmonic differential on the dual mesh for the given
	 * cycle, i.e. the harmonic differential, whose integral along the cycle is
	 * 1.
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
	> DoubleMatrix1D getStandardHarmonicFormOnDualMesh(
				HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
				List<E> cycle) {

		Map<Integer, Integer> tau = getIdentificationMapOnDualMesh(hds, cycle);

		// Consider the surface hds being obtained by identifying two boundary
		// cycles of another surface M via the map above. On this surface we can
		// get a harmonic function with boundary conditions 0 and 1 constant on
		// the identified cycles.
		DoubleMatrix1D h = getStandardHarmonicFunctionOnDualMesh(hds, adapters, cycle,
				tau);

		// The differential of h can be defined on hds.
		DoubleMatrix1D dh = DoubleFactory1D.dense.make(hds.numEdges() / 2);
		Set<F> boundaryVertexSet = EdgeUtility.getDualVertexSet(cycle);

		// build the differences for each edge
		for (E e : hds.getPositiveEdges()) {

			int k = adapters.get(EdgeIndex.class, e, Integer.class);
			EdgeStatus status = EdgeUtility.getDualEdgeStatus(e, cycle,boundaryVertexSet);

			switch (status) {
			case startsAtRightCycle:
				dh.set(k,
						h.get(e.getLeftFace().getIndex())
								- h.get(tau.get(e.getRightFace().getIndex())));
				break;
			case endsAtRightCycle:
				dh.set(k,
						h.get(tau.get(e.getLeftFace().getIndex()))
								- h.get(e.getRightFace().getIndex()));
				break;
			default:
				dh.set(k,
						h.get(e.getLeftFace().getIndex())
								- h.get(e.getRightFace().getIndex()));
				break;
			}
		}

		return dh;
	}

	/**
	 * Returns a harmonic function, which is two valued on the cycle (constant 0
	 * on the left and 1 on the right side).
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
	> DoubleMatrix1D getStandardHarmonicFunctionOnPrimalMesh(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle, Map<Integer, Integer> tau) {
	
		double[] bc1 = new double[tau.size()];
		double[] bc2 = new double[tau.size()];
		for (int i = 0; i < tau.size(); i++) {
			bc2[i] = 1.;
		}
		return getPrimalHarmonicFunction(hds, adapters, cycle, tau, bc1, bc2);
		
	}

	/**
	 * Returns a harmonic function on the dual surface, which is two valued on
	 * the cycle (constant 0 on the left and 1 on the right side).
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
	> DoubleMatrix1D getStandardHarmonicFunctionOnDualMesh(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle, Map<Integer, Integer> tau) {
	
		double[] bc1 = new double[tau.size()];
		double[] bc2 = new double[tau.size()];
		for (int i = 0; i < tau.size(); i++) {
			bc2[i] = 1.;
		}
		return getDualHarmonicFunction(hds, adapters, cycle, tau, bc1, bc2);
		
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
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix1D getPrimalHarmonicFunction(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle, Map<Integer, Integer> tau,
			double[] boundaryCondition1, double[] boundaryCondition2) {

		Set<V> vertexSet = EdgeUtility.getPrimalVertexSet(cycle);
		
		int n = hds.numVertices() + vertexSet.size();
		DoubleMatrix2D laplaceop = LaplaceUtility.getPrimalLaplacian(hds, adapters, cycle, tau);

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

		// finally comment out
		// SimpleMatrixPrintUtility.print(laplaceop,2);
		// System.err.println();
		// System.err.println("determinant: "+dalgebra.det(laplaceop));
		// System.err.println();
		// SimpleMatrixPrintUtility.print(b,2);

		DoubleIterativeSolver solver;
		// solver = new DoubleBiCG(x);
		// solver = new DoubleBiCGstab(x);
		// solver= new DoubleCGLS();
		solver = new DoubleGMRES(x);

		DefaultDoubleIterationMonitor monitor = new DefaultDoubleIterationMonitor();

		// DoubleIterationMonitor monitor= new CGLSDoubleIterationMonitor();
		// configure monitor
		monitor.setMaxIterations(maxIterations);
		monitor.setAbsoluteTolerance(eps);
		monitor.setRelativeTolerance(eps);

		monitor.setDivergenceTolerance(1);

		monitor.setNormType(Norm.Two);

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

		DoubleMatrix1D H = DoubleFactory1D.dense.make(n);

		for (int i = 0; i < n; i++) {
			H.set(i, x.get(i));
		}

		// for (Integer I : tau.keySet()) {
		// System.err.println("h+ = " + H.get(I) + " ; h- = "
		// + H.get(tau.get(I)));
		// }

		return H;
	}

	/**
	 * Returns a harmonic function on the dual surface cut along the given cycle. The
	 * cycle splits in a left and a right side, where boundary conditions can be
	 * specified. Tau is the identification map to obtain the surface given by
	 * the dual hds.
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
	private static <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> DoubleMatrix1D getDualHarmonicFunction(
			HalfEdgeDataStructure<V, E, F> hds, AdapterSet adapters,
			List<E> cycle, Map<Integer, Integer> tau,
			double[] boundaryCondition1, double[] boundaryCondition2) {

		Set<F> vertexSet = EdgeUtility.getDualVertexSet(cycle);
		int n = hds.numFaces() + vertexSet.size();
		
		DoubleMatrix2D laplaceop = LaplaceUtility.getDualLaplacian(hds, adapters, cycle, tau);

		DoubleMatrix1D x = DoubleFactory1D.dense.make(n);

		double[] bcond = new double[n];
		int bcIndex = 0;
		for (F f : vertexSet) {
			int id = f.getIndex();
			bcond[id] = boundaryCondition1[bcIndex];
			bcond[tau.get(id)] = boundaryCondition2[bcIndex];
			bcIndex++;
		}

		DoubleMatrix1D b = DoubleFactory1D.dense.make(bcond);

		// // finally comment out
		// SimpleMatrixPrintUtility.print(laplaceop,2);
		// System.err.println();
		// System.err.println("determinant: "+dalgebra.det(laplaceop));
		// System.err.println();
		// SimpleMatrixPrintUtility.print(b,2);

		DoubleIterativeSolver solver;
		// solver = new DoubleBiCG(x);
		// solver = new DoubleBiCGstab(x);
		// solver= new DoubleCGLS();
		solver = new DoubleGMRES(x);

		DefaultDoubleIterationMonitor monitor = new DefaultDoubleIterationMonitor();

		// DoubleIterationMonitor monitor= new CGLSDoubleIterationMonitor();
		 
		// configure monitor
		monitor.setMaxIterations(maxIterations);
		monitor.setAbsoluteTolerance(eps);
		monitor.setRelativeTolerance(eps);

		monitor.setDivergenceTolerance(1);

		monitor.setNormType(Norm.Two);

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

		DoubleMatrix1D H = DoubleFactory1D.dense.make(n);

		for (int i = 0; i < n; i++) {
			H.set(i, x.get(i));
		}

		// for (Integer I : tau.keySet()) {
		// System.err.println("h+ = " + H.get(I) + " ; h- = "
		// + H.get(tau.get(I)));
		// }

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
	> Map<Integer, Integer> getIdentificationMapOnPrimalMesh(
				HalfEdgeDataStructure<V, E, F> hds, List<E> cycle) {
			
		Map<Integer, Integer> map = new HashMap<Integer, Integer>(100);

		// Get the vertices contained in the cycle.
		Set<V> vertexSet = EdgeUtility.getPrimalVertexSet(cycle);
		int n = hds.numVertices();

		// Build up the map.
		for (V v:vertexSet) {
			map.put(v.getIndex(), n++);
		}
		return map;
	}

	/**
	 * Returns a map which maps the indices of the vertices v0,v1,... contained
	 * in the cycle to the integers n,n+1,... , where n is the number of
	 * vertices of the dual surface (i.e. faces).
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
	> Map<Integer, Integer> getIdentificationMapOnDualMesh(
				HalfEdgeDataStructure<V, E, F> hds, List<E> cycle) {
			
		Map<Integer, Integer> map = new HashMap<Integer, Integer>(100);

		// Get the vertices contained in the cycle.
		Set<F> vertexSet = EdgeUtility.getDualVertexSet(cycle);
		// for (F f : vertexSet) {
		// System.err.print(f.getIndex() + ", ");
		// }
		// System.err.println();
		int n = hds.numFaces();

		// Build up the map.
		for (F f:vertexSet) {
			map.put(f.getIndex(), n++);
		}
		return map;
	}

}
