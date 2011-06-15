package de.varylab.discreteconformal.unwrapper.numerics;

import static de.varylab.jtao.TaoAppAddHess.PreconditionerType.SAME_NONZERO_PATTERN;
import de.varylab.discreteconformal.functional.HyperbolicCircularHolesFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.jpetsc.Mat;
import de.varylab.jpetsc.Vec;
import de.varylab.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.varylab.jtao.TaoAppAddHess;
import de.varylab.jtao.TaoApplication;

public class CHyperbolicApplication extends TaoApplication implements
		TaoAppAddCombinedObjectiveAndGrad, TaoAppAddHess {

	private CoHDS
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
	private HyperbolicCircularHolesFunctional<CoVertex, CoEdge, CoFace> 
		functional = new HyperbolicCircularHolesFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);

	public CHyperbolicApplication(CoHDS hds) {
		this.hds = hds;
	}
	
	@Override
	public double evaluateObjectiveAndGradient(Vec x, Vec g) {
		TaoDomain u = new TaoDomain(x);
		ConformalEnergy E = new ConformalEnergy();
		TaoGradient G = new TaoGradient(g);
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
		return SAME_NONZERO_PATTERN;
	}

	
	public int getDomainDimension() {
		return functional.getDimension(hds);
	}
	
	
}
