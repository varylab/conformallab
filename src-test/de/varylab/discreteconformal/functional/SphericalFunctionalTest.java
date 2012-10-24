package de.varylab.discreteconformal.functional;

import java.util.Random;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.jreality.math.Pn;
import de.jreality.plugin.JRViewer;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.jtem.halfedgetools.functional.MyGradient;
import de.jtem.halfedgetools.functional.MyHessian;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.util.TestUtility;
import de.varylab.discreteconformal.util.UnwrapUtility;

public class SphericalFunctionalTest extends FunctionalTest<CoVertex, CoEdge, CoFace> {

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
	private SphericalFunctional<CoVertex, CoEdge, CoFace>
		functional = new SphericalFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);
	private Random 
		rnd = new Random(); 
	
	
	public static void main(String[] args) {
		JRViewer v = new JRViewer();
		v.addContentUI();
		v.addBasicUI();
		HalfedgeInterface hif = new HalfedgeInterface();
		v.registerPlugin(hif);
		v.startup();
		
		CoHDS hds = new CoHDS(); 
		AdapterSet aSet = new ConformalAdapterSet();
		createOctahedron(hds, aSet);
		hif.addAdapter(new CoPositionAdapter(), true);
		hif.set(hds);
	}
	
	@Override
	public void init() {
		rnd.setSeed(1);
		
		CoHDS hds = new CoHDS(); 
		AdapterSet aSet = new ConformalAdapterSet();
		createOctahedron(hds, aSet);
		
		Random rnd = new Random();
		rnd.setSeed(1);
		for (CoVertex v : hds.getVertices()) {
			Pn.setToLength(v.P, v.P, 0.5 + 1E-2*rnd.nextDouble(), Pn.EUCLIDEAN);
		}
		int n = hds.numVertices();

		// we need small length at start for corresponding spherical length to exist
		Vector x = new DenseVector(n);
		for (int i = 0; i < n; i++) x.set(i, 0);
		MyDomainValue u = new MyDomainValue(x);
		
		UnwrapUtility.prepareInvariantDataHyperbolicAndSpherical(functional, hds, aSet, u);
		
		setFunctional(functional);
		setHDS(hds);
		setXGradient(u);
		setXHessian(u);
		setEps(eps);
		setError(error);
		
		MyGradient g = new MyGradient(new DenseVector(n));
		functional.evaluate(hds, u, null, g, null);
		System.out.println("hand coded: \n" + g);
		
		MyGradient gfd = new MyGradient(new DenseVector(n));
		TestUtility.calculateFDGradient(hds, functional, n, u, gfd);
		System.out.println("finite differences: \n" + gfd);
		
		MyHessian H = new MyHessian(new DenseMatrix(n, n));
		functional.evaluate(hds, u, null, null, H);
		System.out.println("hand coded: \n" + H);
		
		MyHessian Hfd = new MyHessian(new DenseMatrix(n, n));
		TestUtility.calculateFDHessian(hds, functional, n, u, Hfd);
		System.out.println("finite differences: \n" + Hfd);
	}
	
}
