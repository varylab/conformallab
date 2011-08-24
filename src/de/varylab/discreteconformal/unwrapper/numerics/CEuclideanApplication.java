package de.varylab.discreteconformal.unwrapper.numerics;

import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.jtem.jtao.TaoAppAddHess;
import de.jtem.jtao.TaoApplication;
import de.varylab.discreteconformal.functional.ConformalFunctional;
import de.varylab.discreteconformal.functional.EuclideanNewFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;

public class CEuclideanApplication extends TaoApplication implements
		TaoAppAddCombinedObjectiveAndGrad, TaoAppAddHess {

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
//	private EuclideanFunctional<CoVertex, CoEdge, CoFace>
//		functional = new EuclideanFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);
	private EuclideanNewFunctional<CoVertex, CoEdge, CoFace>
		functional = new EuclideanNewFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);
		

	public CEuclideanApplication(CoHDS hds) {
		this.hds = hds;
	}
	
	public ConformalFunctional<CoVertex, CoEdge, CoFace> getFunctional() {
		return functional;
	}

	@Override
	public double evaluateObjectiveAndGradient(Vec x, Vec g) {
		TaoDomain u = new TaoDomain(x);
		TaoGradient G = new TaoGradient(g);
		ConformalEnergy E = new ConformalEnergy();
		functional.evaluate(hds, u, E, G, null);
		g.assemble();
		return E.get();
	}

	@Override
	public PreconditionerType evaluateHessian(Vec x, Mat H, Mat Hpre) {
		TaoDomain u = new TaoDomain(x);
		TaoHessian taoHess = new TaoHessian(H);
		functional.evaluate(hds, u, null, null, taoHess);
		H.assemble();
		return PreconditionerType.SAME_NONZERO_PATTERN;
	}

	
	public int getDomainDimension() {
		return functional.getDimension(hds);
	}
	
	
}
