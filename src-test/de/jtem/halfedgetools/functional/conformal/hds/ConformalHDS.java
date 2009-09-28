package de.jtem.halfedgetools.functional.conformal.hds;

import static java.lang.Math.PI;
import static java.lang.Math.log;

import java.util.Collection;
import java.util.HashSet;

import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.halfedgetools.functional.MyEnergy;
import de.jtem.halfedgetools.functional.conformal.ConformalEuclideanFunctional;
import de.jtem.halfedgetools.functional.conformal.ConformalHyperbolicFunctional;
import de.jtem.halfedgetools.functional.conformal.ConformalAdapters.InitialEnergy;
import de.jtem.halfedgetools.functional.conformal.hds.MyConformalAdapters.CAlpha;
import de.jtem.halfedgetools.functional.conformal.hds.MyConformalAdapters.CLambda;
import de.jtem.halfedgetools.functional.conformal.hds.MyConformalAdapters.CTheta;
import de.jtem.halfedgetools.functional.conformal.hds.MyConformalAdapters.CVariable;

public class ConformalHDS extends HalfEdgeDataStructure<MyConformalVertex, MyConformalEdge, MyConformalFace> {

	private boolean
		texCoordinatesValid = false;
	
	public ConformalHDS() {
		super(MyConformalVertex.class, MyConformalEdge.class, MyConformalFace.class);
	}

	
	/**
	 * Compute algorithm invariant data. Boundary is the natural mesh boundary.
	 * @return the dimension of the parameter space
	 */
	public int prepareInvariantDataEuclidean() {
		HashSet<MyConformalVertex> b = new HashSet<MyConformalVertex>();
		b.addAll(HalfEdgeUtils.boundaryVertices(this));
		return prepareInvariantDataEuclidean(b);
	}
	
	
	/**
	 * Compute algorithm invariant data. Boundary is the natural mesh boundary.
	 * @return the dimension of the parameter space
	 */
	public int prepareInvariantDataHyperbolic() {
		HashSet<MyConformalVertex> b = new HashSet<MyConformalVertex>();
		b.addAll(HalfEdgeUtils.boundaryVertices(this));
		return prepareInvariantDataHyperbolic(b);
	}

	
	/**
	 * Compute algorithm invariant data
	 * @param boundary the boundary vertices which do not belong to the solver system
	 * @return the dimension of the parameter space
	 */
	public int prepareInvariantDataEuclidean(Collection<MyConformalVertex> boundary) {
		// set initial lambdas
		for (final MyConformalEdge e : getPositiveEdges()) {
			final double l = e.getLength();
			e.setLambda(log(l));
			e.getOppositeEdge().setLambda(e.getLambda());
		}
		// set thetas and solver indices
		int dim = 0;
		for (final MyConformalVertex v : getVertices()) {
			if (boundary.contains(v)) {
				v.setTheta(0.0);
				v.setSolverIndex(-1);
			} else {
				v.setTheta(2 * PI);
				v.setSolverIndex(dim++);
			}
		}
		// initial Euclidean energy
		ZeroU zeroU = new ZeroU();
		CVariable var = new CVariable();
		CLambda lambda = new CLambda();
		CAlpha alpha = new CAlpha();
		CTheta theta = new CTheta();
		ZeroInitialEnergy zeroEnergy = new ZeroInitialEnergy();
		MyEnergy E = new MyEnergy();
		ConformalEuclideanFunctional<MyConformalVertex, MyConformalEdge, MyConformalFace, ZeroU> func = new ConformalEuclideanFunctional<MyConformalVertex, MyConformalEdge, MyConformalFace, ZeroU>(var, theta, lambda, alpha, zeroEnergy);
		for (final MyConformalFace f : getFaces()) {
			E.setZero();
			func.triangleEnergyAndAlphas(this, zeroU, f, E);
			f.setInitialEnergy(E.get());
		}
		return dim;
	}
	
	/**
	 * Compute algorithm invariant data
	 * @param boundary the boundary vertices which do not belong to the solver system
	 * @return the dimension of the parameter space
	 */
	public int prepareInvariantDataHyperbolic(Collection<MyConformalVertex> boundary) {
		// set initial lambdas
		for (final MyConformalEdge e : getPositiveEdges()) {
			final double l = e.getLength();
			e.setLambda(2 * log(l));
			e.getOppositeEdge().setLambda(e.getLambda());
		}
		// set thetas and solver indices
		int dim = 0;
		for (final MyConformalVertex v : getVertices()) {
			if (boundary.contains(v)) {
				v.setTheta(0.0);
				v.setSolverIndex(-1);
			} else {
				v.setTheta(2 * PI);
				v.setSolverIndex(dim++);
			}
		}
		// initial hyperbolic energy
		ZeroU zeroU = new ZeroU();
		CVariable var = new CVariable();
		CLambda lambda = new CLambda();
		CTheta theta = new CTheta();
		CAlpha alpha = new CAlpha();
		ZeroInitialEnergy zeroEnergy = new ZeroInitialEnergy();
		MyEnergy E = new MyEnergy();
		ConformalHyperbolicFunctional<MyConformalVertex, MyConformalEdge, MyConformalFace, ZeroU> func = new ConformalHyperbolicFunctional<MyConformalVertex, MyConformalEdge, MyConformalFace, ZeroU>(var, theta, lambda, alpha, zeroEnergy);
		for (final MyConformalFace f : getFaces()) {
			E.setZero();
			func.triangleEnergyAndAlphas(zeroU, f, E);
			f.setInitialEnergy(E.get());
		}
		return dim;
	}
	
	
	
	private class ZeroInitialEnergy implements InitialEnergy<MyConformalFace> {
		@Override
		public double getInitialEnergy(MyConformalFace f) {
			return 0.0;
		}
	}
	
	private class ZeroU implements DomainValue {
		@Override
		public void add(int i, double value) {
		}
		@Override
		public void set(int i, double value) {
		}
		@Override
		public void setZero() {
		}
		@Override
		public double get(int i) {
			return 0.0;
		}
	}
	
	
	
	public boolean isTexCoordinatesValid() {
		return texCoordinatesValid;
	}

	public void setTexCoordinatesValid(boolean texCoordinatesValid) {
		this.texCoordinatesValid = texCoordinatesValid;
	}
	
}
