package de.varylab.discreteconformal.unwrapper.numerics;

import static de.varylab.jtao.TaoAppAddHess.PreconditionerType.SAME_NONZERO_PATTERN;
import de.varylab.discreteconformal.functional.EuclideanCircularHolesFunctional;
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
	private EuclideanCircularHolesFunctional<CoVertex, CoEdge, CoFace>
		functional = new EuclideanCircularHolesFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);
		

	public CEuclideanApplication(CoHDS hds) {
		this.hds = hds;
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
		return SAME_NONZERO_PATTERN;
	}

	
	public int getDomainDimension() {
		// TODO: use method from functional
		int dim = 0;
		for (CoVertex v : hds.getVertices()) {
			if (v.getSolverIndex() >= 0) {
				dim++;
			}
		}
		for (CoEdge e : hds.getPositiveEdges()) {
			if (e.getSolverIndex() >= 0) {
				dim++;
			}
		}
		return dim;
	}
	
	
}
