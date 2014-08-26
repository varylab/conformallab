package de.varylab.discreteconformal.unwrapper.numerics;

import java.util.logging.Logger;

import de.jtem.jpetsc.Vec;
import de.jtem.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.jtem.jtao.TaoApplication;
import de.varylab.discreteconformal.functional.HyperIdealFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CBeta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;

public class CHyperIdealApplication extends TaoApplication implements TaoAppAddCombinedObjectiveAndGrad {

	private Logger
		log = Logger.getLogger(CHyperIdealApplication.class.getName());
	private CoHDS
		hds = null;
	private CTheta
		theta = new CTheta();
	private CVariable
		variable = new CVariable();
	private CAlpha
		alpha = new CAlpha();
	private CBeta
		beta = new CBeta();
	private HyperIdealFunctional<CoVertex, CoEdge, CoFace>
		functional = new HyperIdealFunctional<>(variable, theta, alpha, beta);

	public CHyperIdealApplication(CoHDS hds) {
		this.hds = hds;
	}
	
	@Override
	public double evaluateObjectiveAndGradient(Vec x, Vec g) {
		log.info("eval@" + x);
		TaoDomain u = new TaoDomain(x);
		SimpleEnergy E = new SimpleEnergy();
		TaoGradient G = new TaoGradient(g);
		functional.evaluate(hds, u, E, G, null);
		g.assemble();
		return E.get();
	}
	
	public int getDomainDimension() {
		return functional.getDimension(hds);
	}

}
