package de.varylab.discreteconformal.unwrapper.numerics;

import static de.varylab.discreteconformal.util.SparseUtility.makeNonZeros;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.halfedgetools.functional.Gradient;
import de.jtem.halfedgetools.functional.Hessian;
import de.varylab.discreteconformal.functional.ConformalEuclideanFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.ConformalEnergy;
import de.varylab.mtjoptimization.Optimizable;

public class CEuclideanOptimizable implements Optimizable {

	private CoHDS
		hds = null;
	private CTheta
		theta = new CTheta();
	private CVariable
		variable = new CVariable();
	private CLambda
		lambda = new CLambda();
	private CInitialEnergy
		energy = new CInitialEnergy();
	private CAlpha
		alpha = new CAlpha();
	private ConformalEuclideanFunctional<CoVertex, CoEdge, CoFace>
		functional = new ConformalEuclideanFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);

	public CEuclideanOptimizable(CoHDS hds) {
		this.hds = hds;
	}
	
	private class MTJU implements DomainValue {

		private Vector
			u = null;
		
		public MTJU(Vector u) {
			this.u = u;
		}
		
		@Override
		public void add(int i, double value) {
			u.add(i, value);
		}

		@Override
		public void set(int i, double value) {
			u.set(i, value);
		}

		@Override
		public void setZero() {
			u.zero();
		}
		
		@Override
		public double get(int i) {
			return u.get(i);
		}
		
	}
	
	
	private class MTJGradient implements Gradient {

		private Vector
			G = null;
		
		public MTJGradient(Vector G) {
			this.G = G;
		}
		
		@Override
		public void add(int i, double value) {
			G.add(i, value);
		}

		@Override
		public void set(int i, double value) {
			G.set(i, value);
		}
		
		@Override
		public void setZero() {
			G.zero();
		}
		
	}
	
	
	private class MTJHessian implements Hessian {
		
		private Matrix
			H = null;
		
		public MTJHessian(Matrix H) {
			this.H = H;
		}

		@Override
		public void add(int i, int j, double value) {
			H.add(i, j, value);
		}

		@Override
		public void set(int i, int j, double value) {
			H.set(i, j, value);
		}
		
		@Override
		public void setZero() {
			H.zero();
		}
		
	}
	
	
	@Override
	public Double evaluate(Vector x, Vector gradient, Matrix hessian) {
		MTJU u = new MTJU(x);
		MTJGradient G = new MTJGradient(gradient);
		MTJHessian H = new MTJHessian(hessian);
		ConformalEnergy E = new ConformalEnergy();
		functional.evaluate(hds, u, E, G, H);
		return E.get();
	}

	@Override
	public Double evaluate(Vector x, Vector gradient) {
		MTJU u = new MTJU(x);
		MTJGradient G = new MTJGradient(gradient);
		ConformalEnergy E = new ConformalEnergy();
		functional.evaluate(hds, u, E, G, null);
		return E.get();
	}

	@Override
	public Double evaluate(Vector x, Matrix hessian) {
		MTJU u = new MTJU(x);
		MTJHessian H = new MTJHessian(hessian);
		ConformalEnergy E = new ConformalEnergy();
		functional.evaluate(hds, u, E, null, H);
		return E.get();
	}

	@Override
	public Double evaluate(Vector x) {
		MTJU u = new MTJU(x);
		ConformalEnergy E = new ConformalEnergy();
		functional.evaluate(hds, u, E, null, null);
		return E.get();
	}

	public Integer getDomainDimension() {
		return functional.getDimension(hds);
	}

	@Override
	public Matrix getHessianTemplate() {
		int dim = getDomainDimension();
		return new CompRowMatrix(dim, dim, makeNonZeros(hds));
	}

}
