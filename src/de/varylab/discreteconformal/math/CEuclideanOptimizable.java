package de.varylab.discreteconformal.math;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import de.varylab.discreteconformal.functional.CEuclideanFuctional;
import de.varylab.discreteconformal.functional.CAdapters.Gradient;
import de.varylab.discreteconformal.functional.CAdapters.Hessian;
import de.varylab.discreteconformal.functional.CAdapters.U;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.CVertex;
import de.varylab.discreteconformal.math.Adapters.CAlpha;
import de.varylab.discreteconformal.math.Adapters.CEnergy;
import de.varylab.discreteconformal.math.Adapters.CLambda;
import de.varylab.discreteconformal.math.Adapters.CTheta;
import de.varylab.discreteconformal.math.Adapters.CVariable;
import de.varylab.discreteconformal.math.optimization.Optimizable;

public class CEuclideanOptimizable implements Optimizable {

	private CHDS
		hds = null;
	private CTheta
		theta = new CTheta();
	private CVariable
		variable = new CVariable();
	private CLambda
		lambda = new CLambda();
	private CEnergy
		energy = new CEnergy();
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
		public void addGradient(int i, double value) {
			G.add(i, value);
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
		public void addHessian(int i, int j, double value) {
			H.add(i, j, value);
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
		CEuclideanFuctional.conformalEnergy(hds, u, E, G, H, variable, theta, lambda, alpha, energy);
		return E[0];
	}

	@Override
	public Double evaluate(Vector x, Vector gradient) {
		double[] E = new double[1];
		MTJU u = new MTJU(x);
		MTJGradient G = new MTJGradient(gradient);
		CEuclideanFuctional.conformalEnergy(hds, u, E, G, null, variable, theta, lambda, alpha, energy);
		return E[0];
	}

	@Override
	public Double evaluate(Vector x, Matrix hessian) {
		double[] E = new double[1];
		MTJU u = new MTJU(x);
		MTJHessian H = new MTJHessian(hessian);
		CEuclideanFuctional.conformalEnergy(hds, u, E, null, H, variable, theta, lambda, alpha, energy);
		return E[0];
	}

	@Override
	public Double evaluate(Vector x) {
		double[] E = new double[1];
		MTJU u = new MTJU(x);
		CEuclideanFuctional.conformalEnergy(hds, u, E, null, null, variable, theta, lambda, alpha, energy);
		return E[0];
	}

	@Override
	public Integer getDomainDimension() {
		return hds.getDomainDimension();
	}

}
