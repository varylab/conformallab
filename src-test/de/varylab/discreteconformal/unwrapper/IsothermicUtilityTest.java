package de.varylab.discreteconformal.unwrapper;

import static java.lang.Math.PI;
import static java.lang.Math.log;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.functional.CPEuclideanFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.circlepattern.CirclePatternUtility;
import de.varylab.discreteconformal.unwrapper.isothermic.IsothermicDelaunay;
import de.varylab.discreteconformal.unwrapper.isothermic.IsothermicUtility;
import de.varylab.discreteconformal.unwrapper.isothermic.IsothermicUtility.OppositeAnglesAdapter;
import de.varylab.discreteconformal.unwrapper.numerics.TaoDomain;
import de.varylab.discreteconformal.util.TestUtility;

public class IsothermicUtilityTest extends FunctionalTest<CoVertex, CoEdge, CoFace> {

	static {
		NativePathUtility.set("native");
		Tao.Initialize();
	}
	
	private final double
		EPS = 1E-3,
		EPS2 = 1E-1;
	public CoHDS
		hds = null;
	public CoEdge
		e0 = null;
	public Map<CoEdge, Double> 
		alphaMap = null;
		
	@Override
	@Before
	public void init() {
		hds = new CoHDS();
		CoVertex v0 = hds.addNewVertex();
		CoVertex v1 = hds.addNewVertex();
		CoVertex v2 = hds.addNewVertex();
		CoVertex v3 = hds.addNewVertex();
		CoVertex v4 = hds.addNewVertex();
		
		CoFace f0 = HalfEdgeUtils.constructFaceByVertices(hds, v0, v1, v2).getLeftFace();
		CoFace f1 = HalfEdgeUtils.constructFaceByVertices(hds, v0, v2, v3).getLeftFace();
		CoFace f2 = HalfEdgeUtils.constructFaceByVertices(hds, v0, v3, v4).getLeftFace();
		CoFace f3 = HalfEdgeUtils.constructFaceByVertices(hds, v0, v4, v1).getLeftFace();

		e0 = HalfEdgeUtils.findEdgeBetweenFaces(f1, f0);
		CoEdge e1 = HalfEdgeUtils.findEdgeBetweenFaces(f2, f1);
		CoEdge e2 = HalfEdgeUtils.findEdgeBetweenFaces(f3, f2);
		CoEdge e3 = HalfEdgeUtils.findEdgeBetweenFaces(f0, f3);
		CoEdge e4 = e0.getNextEdge();
		CoEdge e5 = e1.getNextEdge();
		CoEdge e6 = e2.getNextEdge();
		CoEdge e7 = e3.getNextEdge();
		
		
		alphaMap = new HashMap<CoEdge, Double>();
		alphaMap.put(e0,  0.0);
		alphaMap.put(e1,  PI/2);
		alphaMap.put(e2,  0.0);
		alphaMap.put(e3,  PI/2);
		alphaMap.put(e4, -PI/4 + EPS2);
		alphaMap.put(e5,  PI/4);
		alphaMap.put(e6, -PI/4);
		alphaMap.put(e7,  PI/4 - EPS2);
		Set<CoEdge> eSet = new HashSet<CoEdge>(alphaMap.keySet());
		for (CoEdge e : eSet) {
			alphaMap.put(e.getOppositeEdge(), alphaMap.get(e));
		}
		
		Map<CoEdge, Double> betaMap = IsothermicUtility.calculateBetasFromAlphas(hds, alphaMap);
		Map<CoEdge, Double> thetaMap = IsothermicUtility.calculateThetasFromBetas(hds, betaMap);
		Map<CoFace, Double> phiMap = IsothermicUtility.calculatePhisFromBetas(hds, betaMap);
		CPEuclideanFunctional<CoVertex, CoEdge, CoFace>
			fun = new CPEuclideanFunctional<CoVertex, CoEdge, CoFace>(thetaMap, phiMap);
		setFunctional(fun);
		Vec xGrad = new Vec(hds.numFaces());
		xGrad.zeroEntries();
		Vec xHess = new Vec(hds.numFaces());
		xHess.zeroEntries();
		TaoDomain xG = new TaoDomain(xGrad);
		TaoDomain xH = new TaoDomain(xHess);
		setXGradient(xG);
		setXHessian(xH);
		
//		TaoGradient G = new TaoGradient(new Vec(fun.getDimension(hds)));
//		fun.evaluate(hds, xG, null, G, null);
//		System.out.println("Grad: " + G.getVec());
		setHDS(hds);
	}
	
	@Test
	public void testCalculateTriangleAngle() throws Exception {
		double a1 = IsothermicUtility.calculateBeta(PI/2, 0, PI/4);
		Assert.assertEquals(PI/2, a1, 1E-10);
		double a2 = IsothermicUtility.calculateBeta(PI/2, 0, -PI/4);
		Assert.assertEquals(PI/2, a2, 1E-10);
		double a3 = IsothermicUtility.calculateBeta(-PI/2 + EPS, PI/2 - EPS, 0);
		Assert.assertEquals(2*EPS, a3, 1E-10);
		double a4 = IsothermicUtility.calculateBeta(-PI/2 + EPS, PI/2 - EPS, PI/2);
		Assert.assertEquals(PI - 2*EPS, a4, 1E-10); // TODO check fails if EPS=0
		double a5 = IsothermicUtility.calculateBeta(-PI/8, -3*PI/8, PI/2);
		Assert.assertEquals(PI/4, a5, 1E-10);
		double a6 = IsothermicUtility.calculateBeta(PI/8, 3*PI/8, PI/2);
		Assert.assertEquals(PI/4, a6, 1E-10);
		double a7 = IsothermicUtility.calculateBeta(-PI/4, PI/4, PI/2);
		Assert.assertEquals(PI/2, a7, 1E-10);
		double a8 = IsothermicUtility.calculateBeta(-PI/4, PI/4, 0);
		Assert.assertEquals(PI/2, a8, 1E-10);
	}
	
	@Test
	public void testCalculateCirclePatternRadii() throws Exception {
		Map<CoEdge, Double> betaMap = IsothermicUtility.calculateBetasFromAlphas(hds, alphaMap);
		Map<CoEdge, Double> thetaMap = IsothermicUtility.calculateThetasFromBetas(hds, betaMap);
		Map<CoFace, Double> phiMap = IsothermicUtility.calculatePhisFromBetas(hds, betaMap);

		Map<CoFace, Double> rhoMap = CirclePatternUtility.calculateCirclePatternRhos(hds, thetaMap, phiMap);
		Assert.assertEquals(0.0, rhoMap.get(hds.getFace(0)), 1E-10);
		Assert.assertEquals(0.0, rhoMap.get(hds.getFace(1)), 1E-10);
		Assert.assertEquals(-0.110740796199, rhoMap.get(hds.getFace(2)), 1E-11);
		Assert.assertEquals(-0.110740796199, rhoMap.get(hds.getFace(3)), 1E-11);
	}
	
	
	@Test
	public void testCreateDelaunayAngleSystem() throws Exception {
		Map<CoEdge, Double> betaMap = IsothermicUtility.calculateBetasFromAlphas(hds, alphaMap);
		OppositeAnglesAdapter oppositeAnglesAdapter = new OppositeAnglesAdapter(betaMap);
		AdapterSet a = new AdapterSet(oppositeAnglesAdapter);

		IsothermicDelaunay.flip(e0, a);
		
		double g1 = betaMap.get(e0.getNextEdge());
		double g2 = betaMap.get(e0.getPreviousEdge());
		double g3 = betaMap.get(e0.getOppositeEdge().getNextEdge());
		double g4 = betaMap.get(e0.getOppositeEdge().getPreviousEdge());
		Assert.assertEquals(PI/4 + EPS2, g1, 1E-8);
		Assert.assertEquals(PI/4 + EPS2, g2, 1E-8);
		Assert.assertEquals(0.0, g3, 1E-8);
		Assert.assertEquals(0.0, g4, 1E-8);
	}
	
	
	
	@Test
	public void testQuasiConformalFactors() throws Exception {
		double checkFactor = Math.PI;
		CoHDS conformal = TestUtility.readOBJ(getClass(), "triangles.obj");
		AdapterSet a = new ConformalAdapterSet();
		Map<CoEdge, Double> edgeLengthMap = new HashMap<CoEdge, Double>();
		Map<CoEdge, Double> texLengthMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : conformal.getPositiveEdges()) {
			double l = a.get(Length.class, e, Double.class);
			edgeLengthMap.put(e, l);
			edgeLengthMap.put(e.getOppositeEdge(), l);
			texLengthMap.put(e, checkFactor*l);
			texLengthMap.put(e.getOppositeEdge(), checkFactor*l);
		}
		
		Map<CoVertex, Double> uMap = IsothermicUtility.calculateQuasiconformalFactors(conformal, edgeLengthMap, texLengthMap);
		
		for (CoVertex v : conformal.getVertices()) {
			double u1 = uMap.get(v);
			Assert.assertEquals(log(checkFactor), u1, 1E-5);
		}
	}
	
	
}
