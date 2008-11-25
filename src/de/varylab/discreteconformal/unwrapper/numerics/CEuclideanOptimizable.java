package de.varylab.discreteconformal.unwrapper.numerics;

import static de.varylab.discreteconformal.heds.util.SparseUtility.makeNonZeros;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import de.jtem.halfedge.functional.Gradient;
import de.jtem.halfedge.functional.Hessian;
import de.jtem.halfedge.functional.conformal.CEuclideanFuctional;
import de.jtem.halfedge.functional.conformal.CAdapters.U;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.CVertex;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.mtjoptimization.Optimizable;

public class CEuclideanOptimizable implements Optimizable {

	private CHDS
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

	public CEuclideanOptimizable(CHDS hds) {
		this.hds = hds;
	}
	
	private class MTJU implements U<CVertex> {

		private Vector
			u = null;
		
		public MTJU(Vector u) {
			this.u = u;
		}
		
		@Override
		public double getU(CVertex v) {
			if (v.getSolverIndex() >= 0) {
				return u.get(v.getSolverIndex());
			} else {
				return 0;
			}
		}

		@Override
		public void setU(CVertex v, double u) {
			if (v.getSolverIndex() >= 0) {
				this.u.set(v.getSolverIndex(), u);
			}
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
		double[] E = new double[1];
		MTJU u = new MTJU(x);
		MTJGradient G = new MTJGradient(gradient);
		MTJHessian H = new MTJHessian(hessian);
		CEuclideanFuctional.conformalEnergyAndGradient(hds, u, E, G, variable, theta, lambda, alpha, energy);
		CEuclideanFuctional.conformalHessian(hds, u, H, variable, lambda);
		return E[0];
	}

	@Override
	public Double evaluate(Vector x, Vector gradient) {
		double[] E = new double[1];
		MTJU u = new MTJU(x);
		MTJGradient G = new MTJGradient(gradient);
		CEuclideanFuctional.conformalEnergyAndGradient(hds, u, E, G, variable, theta, lambda, alpha, energy);
		return E[0];
	}

	@Override
	public Double evaluate(Vector x, Matrix hessian) {
		double[] E = new double[1];
		MTJU u = new MTJU(x);
		MTJHessian H = new MTJHessian(hessian);
		CEuclideanFuctional.conformalHessian(hds, u, H, variable, lambda);
		return E[0];
	}

	@Override
	public Double evaluate(Vector x) {
		double[] E = new double[1];
		MTJU u = new MTJU(x);
		CEuclideanFuctional.conformalEnergyAndGradient(hds, u, E, null, variable, theta, lambda, alpha, energy);
		return E[0];
	}

	public Integer getDomainDimension() {
		int dim = 0;
		for (CVertex v : hds.getVertices()) {
			if (v.getSolverIndex() >= 0) {
				dim++;
			}
		}
		return dim;
	}

	@Override
	public Matrix getHessianTemplate() {
		int dim = getDomainDimension();
		return new CompRowMatrix(dim, dim, makeNonZeros(hds));
	}

}
