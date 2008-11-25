package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;
import static java.lang.Math.log;

import java.util.Collection;
import java.util.HashSet;

import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.functional.conformal.CEuclideanFuctional;
import de.jtem.halfedge.functional.conformal.CAdapters.InitialEnergy;
import de.jtem.halfedge.functional.conformal.CAdapters.U;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;

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
	public int prepareInvariantData() {
		HashSet<CVertex> b = new HashSet<CVertex>();
		b.addAll(HalfEdgeUtils.boundaryVertices(this));
		return prepareInvariantData(b);
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
	
	
	/**
	 * Compute algorithm invariant data
	 * @param boundary the boundary vertices which do not belong to the solver system
	 * @return the dimension of the parameter space
	 */
	public int prepareInvariantData(Collection<CVertex> boundary) {
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
		ZeroU zeroU = new ZeroU();
		CVariable var = new CVariable();
		CLambda lambda = new CLambda();
		CAlpha alpha = new CAlpha();
		ZeroInitialEnergy zeroEnergy = new ZeroInitialEnergy();
		for (final CFace f : getFaces()) {
			double E = CEuclideanFuctional.triangleEnergyAndAlphas(this, zeroU, f, var, lambda, alpha, zeroEnergy);
			f.setInitialEnergy(E);
		}
		return dim;
	}
	
	
	public boolean isTexCoordinatesValid() {
		return texCoordinatesValid;
	}

	public void setTexCoordinatesValid(boolean texCoordinatesValid) {
		this.texCoordinatesValid = texCoordinatesValid;
	}
	
}
