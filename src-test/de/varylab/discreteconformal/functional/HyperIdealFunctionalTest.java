package de.varylab.discreteconformal.functional;

import static java.lang.Math.abs;

import java.util.Random;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.jtem.halfedgetools.functional.MyEnergy;
import de.jtem.halfedgetools.functional.MyGradient;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.logging.LoggingUtility;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CBeta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;

public class HyperIdealFunctionalTest extends FunctionalTest<CoVertex, CoEdge, CoFace> {

	public static final Double
		eps = 1E-5,
		error = 1E-4;
	private CTheta
		theta = new CTheta();
	private CVariable
		variable = new CVariable();
	private CAlpha
		alpha = new CAlpha();
	private CBeta
		beta = new CBeta();
	private HyperIdealFunctional<CoVertex, CoEdge, CoFace>
		functional = new HyperIdealFunctional<>(variable, theta, alpha, beta);
	
	@BeforeClass
	public static void beforeClass() {
		LoggingUtility.initLogging();
	}
	
	@Override
	public void init() {
		setEps(eps);
		setError(error);
		setFunctional(functional);
		
		CoHDS hds = HyperIdealGenerator.createLawsonSquareTiled();
		int n = functional.getDimension(hds);
		Random rnd = new Random(); 
		rnd.setSeed(1);
		Vector x = new DenseVector(n);
		for (Integer i = 0; i < x.size(); i++) {
			x.set(i, 0.5 + abs(rnd.nextDouble()));
		}
		MyDomainValue u = new MyDomainValue(x);
		setFunctional(functional);
		setHDS(hds);
		setXGradient(u);
	}
	
	@Test@Ignore@Override
	public void testHessian() throws Exception {
	}

	@Test
	public void testGradientWithHyperIdealAndIdealPoints() throws Exception {
		CoHDS hds = HyperIdealGenerator.createLawsonSquareTiledWithBranchPoints();
		int n = functional.getDimension(hds);
		Random rnd = new Random(); 
		rnd.setSeed(1);
		Vector x = new DenseVector(n);
		for (Integer i = 0; i < x.size(); i++) {
			x.set(i, 0.5 + abs(rnd.nextDouble()));
		}
		MyDomainValue u = new MyDomainValue(x);
		setHDS(hds);
		setXGradient(u);
		super.testGradient();
	}
	
	
	@Test
	public void testGradientInTheExtendedDomain() throws Exception {
		CoHDS hds = HyperIdealGenerator.createLawsonSquareTiledWithBranchPoints();
		int n = functional.getDimension(hds);
		Random rnd = new Random(); 
		rnd.setSeed(1);
		Vector x = new DenseVector(n);
		for (Integer i = 0; i < x.size(); i++) {
			x.set(i, 1.2 + abs(rnd.nextDouble()));
		}
		MyDomainValue u = new MyDomainValue(x);
		setHDS(hds);
		setXGradient(u);
		super.testGradient();
	}
	
	@Test
	public void testGradientWithHyperellipticCurve() throws Exception {
		CoHDS hds = HyperIdealGenerator.createLawsonSquareTiledWithBranchPoints();
		int n = functional.getDimension(hds);
		Random rnd = new Random(); 
		rnd.setSeed(1);
		Vector x = new DenseVector(n);
		for (Integer i = 0; i < x.size(); i++) {
			x.set(i, 1.2 + abs(rnd.nextDouble()));
		}
		MyDomainValue u = new MyDomainValue(x);
		setHDS(hds);
		setXGradient(u);
		super.testGradient();
	}
	
	@Test
	public void testFunctionalAtNaNValue() throws Exception {
		CoHDS hds = HyperIdealGenerator.createLawsonHyperelliptic();
		int n = functional.getDimension(hds);
//		Vector x = new DenseVector(new double[] {
//			35.745636463350856, 34.96148036044892, 35.42198104328744, 35.879749523514946, 36.80900102594767, 37.23656501834582, 
//			-83.28921375677265, -39.192883158550664, -36.66444126042561, -38.672883879888026, -34.853144385096, -83.55885278041129, 
//			-83.01782665305134, -39.87814530989889, -33.010547711215295, -33.88691809354694, -82.96421031564259, -40.68162729530564, 
//			-82.28077664119596, -44.64397095058442, -30.486527461909496, -37.55928034865052, -82.8147943723724, -37.82552232466116, 
//			-43.26783038245712, -30.55191678699538, -82.23273951387722, -37.14852659368913, -82.7477086525954, -38.36860352723525, 
//			-33.68458798774745, -82.87230257969732, -40.18581867851664, -39.70844199023845, -33.56093544751282, -82.86828028766375, 
//			-84.0680395032582, -39.13212529085871, -34.594000030866816, -37.45650631253359, -35.87792808148576, -83.04348157235552, 
//			-36.59966232938285, -81.80507713141549, -38.96912737185601, -34.60368564103787, -82.90077030511445, -38.58018732562155, 
//			-39.42404131677899, -33.618237582162614, -83.02830355598282, -40.481302003431075, -33.6222132226343, -82.88592458089232, 
//			-42.910948183197334, -31.02394067870474, -83.1706906471228, -38.424072269634934, -37.0564634872885, -83.25716183867725, 
//			-82.0322281281163, -45.03288212070683, -29.79068034245983, -37.76624515588441, -81.97260700855927, -38.18743539922481, 
//			-83.78324172900801, -39.13314521882563, -33.98623997785822, -39.908962376835525, -33.947545948855485, -82.72513123624627, 
//			-83.1473700954886, -39.306851475389124, -34.42945196285049, -36.72557880308985, -84.06990335805452, -38.08916585400834
//		});
//		Vector x = new DenseVector(new double[]{
//			19.716471454862624, 16.66746398857239, 16.99155189203318, 17.18201781760467, 19.62485697892269, 3.1389175706811905, 
//			-32.17569174060854, -10.18718068377601, -18.739804315407095, -10.277992702726166, -18.146070512594793, -28.440202205036144, 
//			-29.521161431404142, -10.916994891411463, -17.476103760187996, -17.689990680842897, -29.4379041586865, -11.094047260174936, 
//			-31.87192664307049, -7.616209075199606, -16.44525425879046, -19.122022115822375, -29.543760135513512, -10.205699384459967, 
//			-6.915520380211094, -16.567952573801723, -31.92224326030292, -19.03998742248223, -32.17149308976549, -9.772817469629, 
//			-17.628187838140065, -29.106389854402565, -11.210388900674667, -10.785511925313646, -17.70242874220526, -29.42296284070591, 
//			-29.868139367184238, -10.473914427095561, -18.03782474252971, -12.060161524940105, -18.530718562243404, -29.013617352760964, 
//			-18.76787601640074, -31.658530816757455, -10.176927632286258, -18.02284912084168, -29.616630287612804, -10.332074550375708, 
//			-10.495893902059482, -17.739419588016062, -29.58050377282786, -11.038282129726312, -17.66965116440401, -29.432080126631078, 
//			-6.814024858433259, -16.617264733840955, -32.244151328062586, -12.674069722279288, -19.020710105725247, -28.454015791577213, 
//			-31.73254886239917, -7.731253198527611, -16.188440788540635, -19.197221247386544, -31.77151669397451, -9.824735864143662, 
//			-28.580367895440958, -10.640177169049789, -17.77107253608379, -10.830557970773363, -17.75172849274396, -29.35282729979932, 
//			-29.53918174156561, -10.444761159091453, -17.98951085471529, -18.80748425864955, -28.605704267189537, -11.439146658429038
//		});
		Vector x = new DenseVector(new double[]{
			19.146697888829472, 21.845094135762345, 19.961993187843575, 16.018140308468404, 18.15758772600286, 16.451152658066906, 
			-35.874833508021354, -15.511764550539954, -17.90921080418267, -12.493674793908406, -12.993461193786924, -38.44580642670292, 
			-35.85781342745783, -14.020474882871149, -12.981805204418539, -15.809433012308757, -37.17482417871455, -13.227705478276633, 
			-35.29491244650119, -17.559779853808518, -17.281958656412385, -19.95872651344428, -34.370665760316705, -21.25831305560334, 
			-16.859091158820007, -17.40465697142365, -36.22576313321774, -19.872640185136838, -35.5944788931962, -21.18263986394452, 
			-15.635334072191284, -35.16703474936984, -14.415271818334515, -13.289075227918573, -13.397912176212278, -37.6017885784259, 
			-33.84723702606553, -11.6986048222051, -12.507897687551615, -14.9506278233321, -15.661024853537867, -34.88463004980254, 
			-14.44937814529186, -32.728792059655504, -15.428099797258916, -12.492922065863585, -35.032718175081506, -11.556764945485247, 
			-15.22134693630062, -14.259625024286832, -34.60235326261644, -13.309741279929376, -18.635829372709935, -33.83034489106652, 
			-16.757595637042172, -17.453969131462884, -35.667137131493284, -12.453517198371873, -16.90071702888551, -34.87438919143876, 
			-33.03468938842525, -17.674823977136523, -17.02514518616256, -17.072196242996963, -35.19450249740522, -21.234558258459185, 
			-32.38885077677489, -17.656562322612405, -18.747807975469332, -13.334121273378289, -12.327907535277111, -35.65803764852514, 
			-35.657989109183866, -13.054674289945451, -12.910991745328193, -15.668240802688107, -34.86910469220432, -14.864157985299979
		});
		Assert.assertEquals(n, x.size());
		MyEnergy E = new MyEnergy();
		MyDomainValue u = new MyDomainValue(x);
		MyGradient G = new MyGradient(new DenseVector(n));
		functional.evaluate(hds, u, E, G, null);
		Assert.assertNotEquals(Double.NaN, E.get());
		setHDS(hds);
		setXGradient(u);
		super.testGradient();
	}
	
}
