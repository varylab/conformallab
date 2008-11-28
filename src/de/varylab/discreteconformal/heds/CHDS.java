package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;
import static java.lang.Math.log;

import java.util.Collection;
import java.util.HashSet;

import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.functional.MyEnergy;
import de.jtem.halfedge.functional.conformal.CEuclideanFunctional;
import de.jtem.halfedge.functional.conformal.CHyperbolicFunctional;
import de.jtem.halfedge.functional.conformal.CAdapters.InitialEnergy;
import de.jtem.halfedge.functional.conformal.CAdapters.U;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.ConformalEnergy;

public class CHDS extends HalfEdgeDataStructure<CVertex, CEdge, CFace> {

	private boolean
		texCoordinatesValid = false;
	
	public CHDS() {
		super(CVertex.class, CEdge.class, CFace.class);
	}

	
	/**
	 * Compute algorithm invariant data. Boundary is the natural mesh boundary.
	 * @return the dimension of the parameter space
	 */
	public int prepareInvariantDataEuclidean() {
		HashSet<CVertex> b = new HashSet<CVertex>();
		b.addAll(HalfEdgeUtils.boundaryVertices(this));
		return prepareInvariantDataEuclidean(b);
	}
	
	
	/**
	 * Compute algorithm invariant data. Boundary is the natural mesh boundary.
	 * @return the dimension of the parameter space
	 */
	public int prepareInvariantDataHyperbolic() {
		HashSet<CVertex> b = new HashSet<CVertex>();
		b.addAll(HalfEdgeUtils.boundaryVertices(this));
		return prepareInvariantDataHyperbolic(b);
	}

	
	/**
	 * Compute algorithm invariant data
	 * @param boundary the boundary vertices which do not belong to the solver system
	 * @return the dimension of the parameter space
	 */
	public int prepareInvariantDataEuclidean(Collection<CVertex> boundary) {
		// set initial lambdas
		for (final CEdge e : getPositiveEdges()) {
			final double l = e.getLength();
			e.setLambda(log(l));
			e.getOppositeEdge().setLambda(e.getLambda());
		}
		// set thetas and solver indices
		int dim = 0;
		for (final CVertex v : getVertices()) {
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
		CEuclideanFunctional<CVertex, CEdge, CFace> func = new CEuclideanFunctional<CVertex, CEdge, CFace>(var, theta, lambda, alpha, zeroEnergy);
		for (final CFace f : getFaces()) {
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
	public int prepareInvariantDataHyperbolic(Collection<CVertex> boundary) {
		// set initial lambdas
		for (final CEdge e : getPositiveEdges()) {
			final double l = e.getLength();
			e.setLambda(2*log(l));
			e.getOppositeEdge().setLambda(e.getLambda());
		}
		// set thetas and solver indices
		int dim = 0;
		for (final CVertex v : getVertices()) {
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
		CHyperbolicFunctional<CVertex, CEdge, CFace> func = new CHyperbolicFunctional<CVertex, CEdge, CFace>(var, theta, lambda, alpha, zeroEnergy);
		for (final CFace f : getFaces()) {
			E.setZero();
			func.triangleEnergyAndAlphas(zeroU, f, E);
			f.setInitialEnergy(E.get());
		}
		return dim;
	}
	
	
	
	private class ZeroInitialEnergy implements InitialEnergy<CFace> {
		@Override
		public double getInitialEnergy(CFace f) {
			return 0.0;
		}
	}
	
	private class ZeroU implements U<CVertex> {
		@Override
		public double getU(CVertex v) {
			return 0;
		}
		@Override
		public void setU(CVertex v, double u) {
		}
	}
	
	
	
	public boolean isTexCoordinatesValid() {
		return texCoordinatesValid;
	}

	public void setTexCoordinatesValid(boolean texCoordinatesValid) {
		this.texCoordinatesValid = texCoordinatesValid;
	}
	
}
