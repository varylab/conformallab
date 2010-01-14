package de.jtem.halfedgetools.functional.conformal;

import java.util.Random;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.jtem.halfedgetools.functional.conformal.hds.ConformalHDS;
import de.jtem.halfedgetools.functional.conformal.hds.MyConformalEdge;
import de.jtem.halfedgetools.functional.conformal.hds.MyConformalFace;
import de.jtem.halfedgetools.functional.conformal.hds.MyConformalVertex;
import de.jtem.halfedgetools.functional.conformal.hds.MyConformalAdapters.CAlpha;
import de.jtem.halfedgetools.functional.conformal.hds.MyConformalAdapters.CInitialEnergy;
import de.jtem.halfedgetools.functional.conformal.hds.MyConformalAdapters.CLambda;
import de.jtem.halfedgetools.functional.conformal.hds.MyConformalAdapters.CTheta;
import de.jtem.halfedgetools.functional.conformal.hds.MyConformalAdapters.CVariable;
import de.varylab.discreteconformal.functional.ConformalHyperbolicFunctional;

public class ConformalHyperbolicFunctionalTest extends FunctionalTest<MyConformalVertex, MyConformalEdge, MyConformalFace> {

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
	private ConformalHyperbolicFunctional<MyConformalVertex, MyConformalEdge, MyConformalFace>
		functional = new ConformalHyperbolicFunctional<MyConformalVertex, MyConformalEdge, MyConformalFace>(variable, theta, lambda, alpha, energy);
	
	
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
		
		setFuctional(functional);
		setHDS(hds);
		setXGradient(u);
		setXHessian(u);
		setEps(eps);
		setError(error);
	}
	
	
}
