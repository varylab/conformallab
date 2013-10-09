package de.varylab.discreteconformal.functional;

import static java.lang.Math.log;

import java.util.Random;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Vector.Norm;

import org.junit.Assert;
import org.junit.Test;

import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
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
import de.varylab.discreteconformal.unwrapper.numerics.MTJGradient;
import de.varylab.discreteconformal.util.UnwrapUtility;
import de.varylab.discreteconformal.util.UnwrapUtility.ZeroU;

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
	
	
	@Override
	public void init() {
		rnd.setSeed(1);
		
		CoHDS hds = new CoHDS(); 
		AdapterSet aSet = new ConformalAdapterSet();
		createOctahedron(hds, aSet);
		
		for (CoVertex v : hds.getVertices()) {
			Pn.setToLength(v.P, v.P, 0.5 + 1E-2*rnd.nextDouble(), Pn.EUCLIDEAN);
		}
		int n = hds.numVertices();

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
//		MyGradient g = new MyGradient(new DenseVector(n));
//		functional.evaluate(hds, u, null, g, null);
//		System.out.println("hand coded: \n" + g);
//		
//		MyGradient gfd = new MyGradient(new DenseVector(n));
//		TestUtility.calculateFDGradient(hds, functional, n, u, gfd);
//		System.out.println("finite differences: \n" + gfd);
//		
//		MyHessian H = new MyHessian(new DenseMatrix(n, n));
//		functional.evaluate(hds, u, null, null, H);
//		System.out.println("hand coded: \n" + H);
//		
//		MyHessian Hfd = new MyHessian(new DenseMatrix(n, n));
//		TestUtility.calculateFDHessian(hds, functional, n, u, Hfd);
//		System.out.println("finite differences: \n" + Hfd);
	}
	
	
	@Test
	public void testReducedGradient() throws Exception {
		functional.setReduced(true);
		testGradient();
	}
	
	@Test
	public void testReducedHessian() throws Exception {
		functional.setReduced(true);
		testHessian();
	}
	
	@Override
	@Test
	public void testGradient() throws Exception {
		functional.setReduced(false);
		super.testGradient();
	}
	
	@Override
	@Test
	public void testHessian() throws Exception {
		functional.setReduced(false);
		super.testHessian();
	}
	
	
	@Test
	public void testCriticalPoint() throws Exception {
		rnd.setSeed(2);
		
		CoHDS hds = new CoHDS(); 
		AdapterSet aSet = new ConformalAdapterSet();
		hds.addNewVertices(12);
		for (CoVertex v : hds.getVertices()) {
			double[] pos = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
			Rn.setEuclideanNorm(pos, 1.0, pos);
			aSet.set(Position.class, v, pos);
		}
		ConvexHull.convexHull(hds, aSet);
		
//		ConverterHeds2JR converter = new ConverterHeds2JR();
//		IndexedFaceSet ifs = converter.heds2ifs(hds, aSet);
//		WriterOBJ.write(ifs, new FileOutputStream("critical.obj"));
		
		for (CoEdge e : hds.getPositiveEdges()) {
			double l = aSet.get(Length.class, e, Double.class);
			double lambda = 2 * log(l / 2);
			e.setLambda(lambda);
			e.getOppositeEdge().setLambda(lambda);
		}
		
		
//		CoHDS hds = new CoHDS(); 
//		HalfEdgeUtils.addIcosahedron(hds);
//		double l = 1 / Math.sin(2 * PI / 5);
//		for (CoEdge e : hds.getEdges()) {
//			e.setLambda(2 * log(l / 2));
//		}

		for (CoVertex v : hds.getVertices()) {
			v.setSolverIndex(v.getIndex());
		}
		
		Vector gVec = new DenseVector(hds.numVertices());
		MTJGradient G = new MTJGradient(gVec);
		ZeroU zeroU = new ZeroU();
		functional.conformalEnergyAndGradient(hds, zeroU, null, G);
		
		// check flatness
//		for (CoVertex v : hds.getVertices()) {
//			double a = 0.0;
//			for (CoEdge e : HalfEdgeUtils.incomingEdges(v)) {
//				a += e.getPreviousEdge().getAlpha();
//			}
//			Assert.assertEquals(2*PI, a, 1E-8);
//		}
		
		System.out.println(gVec);
		// check critical point
		Assert.assertEquals(0.0, gVec.norm(Norm.Two), 1E-8);
	}
	
}
