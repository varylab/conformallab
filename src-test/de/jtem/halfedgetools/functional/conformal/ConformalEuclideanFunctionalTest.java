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
import de.varylab.discreteconformal.functional.ConformalEuclideanFunctional;

public class ConformalEuclideanFunctionalTest extends FunctionalTest<MyConformalVertex, MyConformalEdge, MyConformalFace> {

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
	public ConformalEuclideanFunctional<MyConformalVertex, MyConformalEdge, MyConformalFace>
		functional = new ConformalEuclideanFunctional<MyConformalVertex, MyConformalEdge, MyConformalFace>(variable, theta, lambda, alpha, energy);
	
	@Override
	public void init() {
		ConformalHDS hds = new ConformalHDS(); 
		createTetrahedron(hds);
		hds.prepareInvariantDataEuclidean();
		
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
	}
	
}
