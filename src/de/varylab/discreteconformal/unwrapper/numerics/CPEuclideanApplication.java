package de.varylab.discreteconformal.unwrapper.numerics;

import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.PETSc;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.jtem.jtao.TaoAppAddHess;
import de.jtem.jtao.TaoApplication;
import de.varylab.discreteconformal.functional.CPEuclideanFunctional;
import de.varylab.discreteconformal.functional.CPEuclideanFunctional.Phi;
import de.varylab.discreteconformal.functional.CPEuclideanFunctional.Theta;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.SparseUtility;

public class CPEuclideanApplication extends TaoApplication implements
		TaoAppAddCombinedObjectiveAndGrad, TaoAppAddHess {

	private CoHDS
		hds = null;
	private CPEuclideanFunctional<CoVertex, CoEdge, CoFace>
		functional = null;
		

	public CPEuclideanApplication(CoHDS hds, Theta<CoEdge> theta, Phi<CoFace> phi) {
		this.hds = hds;
		this.functional = new CPEuclideanFunctional<CoVertex, CoEdge, CoFace>(theta, phi);
	}
	
	public CPEuclideanFunctional<CoVertex, CoEdge, CoFace> getFunctional() {
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
	
	public Mat getHessianTemplate() {
		int dim = getDomainDimension();
		int[][] sparceStructure = functional.getNonZeroPattern(hds);
		int[] nonZeros = SparseUtility.getPETScNonZeros(sparceStructure);
		Mat H = Mat.createSeqAIJ(dim, dim, PETSc.PETSC_DEFAULT, nonZeros);
		H.assemble();
		return H;
	}
	
	
}
