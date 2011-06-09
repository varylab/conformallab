package de.varylab.discreteconformal.unwrapper.numerics;

import static de.varylab.discreteconformal.util.SparseUtility.makeNonZeros;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import de.varylab.discreteconformal.functional.HyperbolicFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.mtjoptimization.Optimizable;

public class CHyperbolicOptimizable implements Optimizable {

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
	private HyperbolicFunctional<CoVertex, CoEdge, CoFace>
		functional = new HyperbolicFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);

	public CHyperbolicOptimizable(CoHDS hds) {
		this.hds = hds;
	}
	
	@Override
	public Double evaluate(Vector x, Vector gradient, Matrix hessian) {
		MTJDomain u = new MTJDomain(x);
		MTJGradient G = new MTJGradient(gradient);
		MTJHessian H = new MTJHessian(hessian);
		ConformalEnergy E = new ConformalEnergy();
		functional.evaluate(hds, u, E, G, H);
		return E.get();
	}

	@Override
	public Double evaluate(Vector x, Vector gradient) {
		MTJDomain u = new MTJDomain(x);
		MTJGradient G = new MTJGradient(gradient);
		ConformalEnergy E = new ConformalEnergy();
		functional.evaluate(hds, u, E, G, null);
		return E.get();
	}

	@Override
	public Double evaluate(Vector x, Matrix hessian) {
		MTJDomain u = new MTJDomain(x);
		MTJHessian H = new MTJHessian(hessian);
		ConformalEnergy E = new ConformalEnergy();
		functional.evaluate(hds, u, E, null, H);
		return E.get();
	}

	@Override
	public Double evaluate(Vector x) {
		MTJDomain u = new MTJDomain(x);
		ConformalEnergy E = new ConformalEnergy();
		functional.evaluate(hds, u, E, null, null);
		return E.get();
	}

	@Override
	public Integer getDomainDimension() {
		return functional.getDimension(hds);
	}

	@Override
	public Matrix getHessianTemplate() {
		int dim = getDomainDimension();
		return new CompRowMatrix(dim, dim, makeNonZeros(hds));
	}

}
