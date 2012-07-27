package de.varylab.discreteconformal.unwrapper.numerics;

import java.util.Map;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.PETSc;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.jtem.jtao.TaoAppAddHess;
import de.jtem.jtao.TaoApplication;
import de.varylab.discreteconformal.functional.CPEuclideanFunctional;
import de.varylab.discreteconformal.util.SparseUtility;

public class CPEuclideanApplication <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>,
	HDS extends HalfEdgeDataStructure<V, E, F>
> extends TaoApplication implements TaoAppAddCombinedObjectiveAndGrad, TaoAppAddHess {

	private HDS
		hds = null;
	private CPEuclideanFunctional<V, E, F>
		functional = null;
		

	public CPEuclideanApplication(HDS hds, Map<E, Double> thetaMap, Map<F, Double> phiMap) {
		this.hds = hds;
		this.functional = new CPEuclideanFunctional<V, E, F>(thetaMap, phiMap);
	}
	
	public CPEuclideanFunctional<V, E, F> getFunctional() {
		return functional;
	}
	
	@Override
	public double evaluateObjectiveAndGradient(Vec x, Vec g) {
		TaoDomain u = new TaoDomain(x);
		TaoGradient G = new TaoGradient(g);
		SimpleEnergy E = new SimpleEnergy();
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
