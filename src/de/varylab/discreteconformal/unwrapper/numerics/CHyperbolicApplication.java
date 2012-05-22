package de.varylab.discreteconformal.unwrapper.numerics;

import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.PETSc;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.jtem.jtao.TaoAppAddHess;
import de.jtem.jtao.TaoApplication;
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
import de.varylab.discreteconformal.util.SparseUtility;

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
	
	public HyperbolicCircularHolesFunctional<CoVertex, CoEdge, CoFace> getFunctional() {
		return functional;
	}
	
	@Override
	public double evaluateObjectiveAndGradient(Vec x, Vec g) {
		TaoDomain u = new TaoDomain(x);
		SimpleEnergy E = new SimpleEnergy();
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
		return PreconditionerType.SAME_NONZERO_PATTERN;
	}

	
	public int getDomainDimension() {
		return functional.getDimension(hds);
	}
	
	public Mat getHessianTemplate() {
		int dim = getDomainDimension();
		int[][] sparceStructure = functional.getNonZeroPattern(hds);
		int[] nonZeros = SparseUtility.getPETScNonZeros(sparceStructure);
		Mat H = Mat.createSeqAIJ(dim, dim, PETSc.PETSC_DEFAULT, nonZeros);
		H.assemble();
		return H;
	}
	
}
