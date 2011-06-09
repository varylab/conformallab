package de.varylab.discreteconformal.functional;

import java.util.Random;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.varylab.discreteconformal.functional.HyperbolicFunctional;
import de.varylab.discreteconformal.functional.hds.ConformalHDS;
import de.varylab.discreteconformal.functional.hds.MyConformalEdge;
import de.varylab.discreteconformal.functional.hds.MyConformalFace;
import de.varylab.discreteconformal.functional.hds.MyConformalVertex;
import de.varylab.discreteconformal.functional.hds.MyConformalAdapters.CAlpha;
import de.varylab.discreteconformal.functional.hds.MyConformalAdapters.CInitialEnergy;
import de.varylab.discreteconformal.functional.hds.MyConformalAdapters.CLambda;
import de.varylab.discreteconformal.functional.hds.MyConformalAdapters.CTheta;
import de.varylab.discreteconformal.functional.hds.MyConformalAdapters.CVariable;

public class HyperbolicFunctionalTest extends FunctionalTest<MyConformalVertex, MyConformalEdge, MyConformalFace> {

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
	private HyperbolicFunctional<MyConformalVertex, MyConformalEdge, MyConformalFace>
		functional = new HyperbolicFunctional<MyConformalVertex, MyConformalEdge, MyConformalFace>(variable, theta, lambda, alpha, energy);
	
	
	@Override
	public void init() {
		ConformalHDS hds = new ConformalHDS(); 
		createCube(hds);
		hds.removeFace(hds.getFace(0));
		hds.prepareInvariantDataHyperbolic();
		
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
