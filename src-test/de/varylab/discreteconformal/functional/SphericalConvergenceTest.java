package de.varylab.discreteconformal.functional;

import static de.varylab.discreteconformal.util.SparseUtility.makeNonZeros;

import java.util.Random;

import junit.framework.Assert;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

import org.junit.Test;

import de.jreality.math.Pn;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.unwrapper.numerics.CSphericalOptimizable;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.mtjoptimization.NotConvergentException;
import de.varylab.mtjoptimization.newton.NewtonOptimizer;
import de.varylab.mtjoptimization.newton.NewtonOptimizer.Solver;

public class SphericalConvergenceTest  {

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
	private SphericalFunctional<CoVertex, CoEdge, CoFace>
		functional = new SphericalFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);
	private Random 
		rnd = new Random();
	
	@Test
	public void testSphericalConvergence() {
		rnd.setSeed(1);
		
		CoHDS hds = new CoHDS(); 
		AdapterSet aSet = new ConformalAdapterSet();
		FunctionalTest.createOctahedron(hds, aSet);
		
		for (CoVertex v : hds.getVertices()) {
			Pn.setToLength(v.P, v.P, 0.5 + 1E-2*rnd.nextDouble(), Pn.EUCLIDEAN);
		}
		
		int n = hds.numVertices();
		Vector x = new DenseVector(n);
		for (int i = 0; i < n; i++) x.set(i, 0);
		MyDomainValue u = new MyDomainValue(x);
		UnwrapUtility.prepareInvariantDataHyperbolicAndSpherical(functional, hds, aSet, u);
		
////		 one edge is circular
//		for (CoEdge e : hds.getPositiveEdges()) {
//			CoVertex s = e.getStartVertex();
//			CoVertex t = e.getTargetVertex();
//			if (isBoundaryVertex(s) || isBoundaryVertex(t)) {
//				continue;
//			}
//			e.info = new CustomEdgeInfo();
//			e.info.circularHoleEdge = true;
//			break;
//		}
		
		CSphericalOptimizable opt = new CSphericalOptimizable(hds);
		// optimization
		DenseVector ux = new DenseVector(n);
		Matrix H = new CompRowMatrix(n,n, makeNonZeros(hds));
		NewtonOptimizer optimizer = new NewtonOptimizer(H);
//	cannot use a step controller here, it stops the gradient flow		
//		optimizer.setStepController(new ArmijoStepController());
		optimizer.setSolver(Solver.CGS);
		optimizer.setError(1E-10);
		optimizer.setMaxIterations(10);
		try {
			optimizer.minimize(ux, opt);
		} catch (NotConvergentException e) {
			Assert.fail(e.getLocalizedMessage());
		}
	}
	
	
}