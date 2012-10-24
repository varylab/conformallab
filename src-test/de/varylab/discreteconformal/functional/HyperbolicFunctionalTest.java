package de.varylab.discreteconformal.functional;

import java.util.Random;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
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
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.UnwrapUtility.ZeroU;

public class HyperbolicFunctionalTest extends FunctionalTest<CoVertex, CoEdge, CoFace> {

	public static final Double
		eps = 1E-5,
		error = 1E-4;
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
	private HyperbolicFunctional<CoVertex, CoEdge, CoFace>
		functional = new HyperbolicFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);
	
	
	@Override
	public void init() {
		CoHDS hds = new CoHDS(); 
		AdapterSet aSet = new ConformalAdapterSet();
		createTriangulatedCube(hds, aSet);
		
		ZeroU zeroU = new ZeroU();
		UnwrapUtility.prepareInvariantDataHyperbolicAndSpherical(functional, hds, aSet, zeroU);
		
		int n = functional.getDimension(hds);
		Random rnd = new Random(); 
		rnd.setSeed(1);
		
		Vector x = new DenseVector(n);
		for (Integer i = 0; i < x.size(); i++) {
			x.set(i, rnd.nextDouble() - 0.5);
		}
		MyDomainValue u = new MyDomainValue(x);
		
		setFunctional(functional);
		setHDS(hds);
		setXGradient(u);
		setXHessian(u);
		setEps(eps);
		setError(error);
	}
	
	
}
