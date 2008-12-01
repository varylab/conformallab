package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;
import static java.lang.Math.log;

import java.util.Collection;
import java.util.HashSet;

import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.functional.DomainValue;
import de.jtem.halfedge.functional.MyEnergy;
import de.jtem.halfedge.functional.conformal.ConformalEuclideanFunctional;
import de.jtem.halfedge.functional.conformal.ConformalHyperbolicFunctional;
import de.jtem.halfedge.functional.conformal.ConformalAdapters.InitialEnergy;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.ConformalEnergy;

public class CoHDS extends HalfEdgeDataStructure<CoVertex, CoEdge, CoFace> {

	private boolean
		texCoordinatesValid = false;
	
	public CoHDS() {
		super(CoVertex.class, CoEdge.class, CoFace.class);
	}

	
	/**
	 * Compute algorithm invariant data. Boundary is the natural mesh boundary.
	 * @return the dimension of the parameter space
	 */
	public int prepareInvariantDataEuclidean() {
		HashSet<CoVertex> b = new HashSet<CoVertex>();
		b.addAll(HalfEdgeUtils.boundaryVertices(this));
		return prepareInvariantDataEuclidean(b);
	}
	
	
	/**
	 * Compute algorithm invariant data. Boundary is the natural mesh boundary.
	 * @return the dimension of the parameter space
	 */
	public int prepareInvariantDataHyperbolic() {
		HashSet<CoVertex> b = new HashSet<CoVertex>();
		b.addAll(HalfEdgeUtils.boundaryVertices(this));
		return prepareInvariantDataHyperbolic(b);
	}

	
	/**
	 * Compute algorithm invariant data
	 * @param boundary the boundary vertices which do not belong to the solver system
	 * @return the dimension of the parameter space
	 */
	public int prepareInvariantDataEuclidean(Collection<CoVertex> boundary) {
		// set initial lambdas
		for (final CoEdge e : getPositiveEdges()) {
			final double l = e.getLength();
			e.setLambda(log(l));
			e.getOppositeEdge().setLambda(e.getLambda());
		}
		// set thetas and solver indices
		int dim = 0;
		for (final CoVertex v : getVertices()) {
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
		ConformalEuclideanFunctional<CoVertex, CoEdge, CoFace, ZeroU> func = new ConformalEuclideanFunctional<CoVertex, CoEdge, CoFace, ZeroU>(var, theta, lambda, alpha, zeroEnergy);
		for (final CoFace f : getFaces()) {
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
	public int prepareInvariantDataHyperbolic(Collection<CoVertex> boundary) {
		// set initial lambdas
		for (final CoEdge e : getPositiveEdges()) {
			final double l = e.getLength();
			e.setLambda(2*log(l));
			e.getOppositeEdge().setLambda(e.getLambda());
		}
		// set thetas and solver indices
		int dim = 0;
		for (final CoVertex v : getVertices()) {
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
		ConformalEnergy E = new ConformalEnergy();
		ConformalHyperbolicFunctional<CoVertex, CoEdge, CoFace, ZeroU> func = new ConformalHyperbolicFunctional<CoVertex, CoEdge, CoFace, ZeroU>(var, theta, lambda, alpha, zeroEnergy);
		for (final CoFace f : getFaces()) {
			E.setZero();
			func.triangleEnergyAndAlphas(zeroU, f, E);
			f.setInitialEnergy(E.get());
		}
		return dim;
	}
	
	
	
	private class ZeroInitialEnergy implements InitialEnergy<CoFace> {
		@Override
		public double getInitialEnergy(CoFace f) {
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
