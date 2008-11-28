package de.varylab.discreteconformal.unwrapper.numerics;

import static de.varylab.jpetsc.Vec.InsertMode.INSERT_VALUES;
import static de.varylab.jtao.TaoAppAddHess.PreconditionerType.SAME_NONZERO_PATTERN;
import de.jtem.halfedge.functional.Gradient;
import de.jtem.halfedge.functional.Hessian;
import de.jtem.halfedge.functional.conformal.CHyperbolicFunctional;
import de.jtem.halfedge.functional.conformal.CAdapters.U;
import de.varylab.discreteconformal.heds.CEdge;
import de.varylab.discreteconformal.heds.CFace;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.CVertex;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.ConformalEnergy;
import de.varylab.jpetsc.Mat;
import de.varylab.jpetsc.Vec;
import de.varylab.jpetsc.Mat.InsertMode;
import de.varylab.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.varylab.jtao.TaoAppAddHess;
import de.varylab.jtao.TaoApplication;

public class CHyperbolicApplication extends TaoApplication implements
		TaoAppAddCombinedObjectiveAndGrad, TaoAppAddHess {

	private CHDS
		hds = null;
	private CTheta
		theta = new CTheta();
	private CVariable
		variable = new CVariable();
	private CLambda
		lambda = new CLambda();
	private CAlpha
		alpha = new CAlpha();
	private CInitialEnergy
		energy = new CInitialEnergy();
	private CHyperbolicFunctional<CVertex, CEdge, CFace> 
		functional = new CHyperbolicFunctional<CVertex, CEdge, CFace>(variable, theta, lambda, alpha, energy);
		

	public CHyperbolicApplication(CHDS hds) {
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
		public void add(int i, double value) {
			G.add(i, value);
		}

		@Override
		public void set(int i, double value) {
			G.setValue(i, value, INSERT_VALUES);
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
		public void add(int i, int j, double value) {
			H.add(i, j, value);
		}

		@Override
		public void setZero() {
			H.zeroEntries();
		}

		@Override
		public void set(int i, int j, double value) {
			H.setValue(i, j, value, InsertMode.INSERT_VALUES);
		}
		
	}
	
	
	@Override
	public double evaluateObjectiveAndGradient(Vec x, Vec g) {
		TaoU u = new TaoU(x);
		ConformalEnergy E = new ConformalEnergy();
		TaoGradient G = new TaoGradient(g);
		functional.evaluate(hds, u, E, G, null);
		g.assemble();
		return E.get();
	}

	@Override
	public PreconditionerType evaluateHessian(Vec x, Mat H, Mat Hpre) {
		TaoU u = new TaoU(x);
		TaoHessian taoHess = new TaoHessian(H);
		functional.evaluate(hds, u, null, null, taoHess);
		H.assemble();
		return SAME_NONZERO_PATTERN;
	}

	
	public int getDomainDimension() {
		return functional.getDimension(hds);
	}
	
	
}
