package de.varylab.discreteconformal.math;

import static de.varylab.jpetsc.Vec.InsertMode.INSERT_VALUES;
import static de.varylab.jtao.TaoAppAddHess.PreconditionerType.SAME_NONZERO_PATTERN;
import de.varylab.discreteconformal.functional.CEuclideanFuctional;
import de.varylab.discreteconformal.functional.CEuclideanFuctional.Gradient;
import de.varylab.discreteconformal.functional.CEuclideanFuctional.Hessian;
import de.varylab.discreteconformal.functional.CEuclideanFuctional.U;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.CVertex;
import de.varylab.discreteconformal.math.Adapters.CAlpha;
import de.varylab.discreteconformal.math.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.math.Adapters.CLambda;
import de.varylab.discreteconformal.math.Adapters.CTheta;
import de.varylab.discreteconformal.math.Adapters.CVariable;
import de.varylab.jpetsc.Mat;
import de.varylab.jpetsc.Vec;
import de.varylab.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.varylab.jtao.TaoAppAddHess;
import de.varylab.jtao.TaoApplication;

public class CEuclideanApplication extends TaoApplication implements
		TaoAppAddCombinedObjectiveAndGrad, TaoAppAddHess {

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
		

	public CEuclideanApplication(CHDS hds) {
		this.hds = hds;
	}


	public static class TaoU implements U<CVertex> {

		private Vec
			u = null;
		
		public TaoU(Vec u) {
			this.u = u;
		}
		
		@Override
		public double getU(CVertex v) {
			if (v.getSolverIndex() >= 0) {
				return u.getValue(v.getSolverIndex());
			} else {
				return 0;
			}
		}

		@Override
		public void setU(CVertex v, double u) {
			if (v.getSolverIndex() >= 0) {
				this.u.setValue(v.getSolverIndex(), u, INSERT_VALUES);
			}
		}
		
	}
	
	
	private static class TaoGradient implements Gradient {

		private Vec
			G = null;
		
		public TaoGradient(Vec G) {
			this.G = G;
		}
		
		@Override
		public void addGradient(int i, double value) {
			G.add(i, value);
		}

		@Override
		public void setZero() {
			G.zeroEntries();
		}
		
	}
	
	
	private static class TaoHessian implements Hessian {
		
		private Mat
			H = null;
		
		public TaoHessian(Mat H) {
			this.H = H;
		}

		@Override
		public void addHessian(int i, int j, double value) {
			H.add(i, j, value);
		}

		@Override
		public void setZero() {
			H.zeroEntries();
		}
		
	}
	

	@Override
	public double evaluateObjectiveAndGradient(Vec x, Vec g) {
		double[] E = new double[1];
		TaoU u = new TaoU(x);
		TaoGradient taoGrad = new TaoGradient(g);
		CEuclideanFuctional.conformalEnergyAndGradient(hds, u, E, taoGrad, variable, theta, lambda, alpha, energy);
		g.assemble();
		return E[0];
	}

	@Override
	public PreconditionerType evaluateHessian(Vec x, Mat H, Mat Hpre) {
		TaoU u = new TaoU(x);
		TaoHessian taoHess = new TaoHessian(H);
		CEuclideanFuctional.conformalHessian(hds, u, taoHess, variable, lambda);
		H.assemble();
		return SAME_NONZERO_PATTERN;
	}

	
	public int getDomainDimension() {
		int dim = 0;
		for (CVertex v : hds.getVertices()) {
			if (v.getSolverIndex() >= 0) {
				dim++;
			}
		}
		return dim;
	}
	
	
}
