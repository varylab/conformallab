package de.varylab.discreteconformal.unwrapper;

import static java.lang.Math.PI;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import cern.colt.Arrays;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.functional.CPEuclideanFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.IsothermicUtility.CPLayoutAdapters;
import de.varylab.discreteconformal.unwrapper.numerics.TaoDomain;
import de.varylab.discreteconformal.unwrapper.numerics.TaoGradient;

public class IsothermicUtilityTest extends FunctionalTest<CoVertex, CoEdge, CoFace> {

	private final double
		EPS = 1E-2;
	public CoHDS
		hds = null;
	public Map<CoEdge, Double> 
		alphaMap = null;
		
	@Before
	public void init() {
		NativePathUtility.set("native");
		Tao.Initialize();
		
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

		CoEdge e0 = HalfEdgeUtils.findEdgeBetweenFaces(f1, f0);
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
		alphaMap.put(e4, -PI/4);
		alphaMap.put(e5,  PI/4);
		alphaMap.put(e6, -PI/4);
		alphaMap.put(e7,  PI/4);
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
		
		TaoGradient G = new TaoGradient(new Vec(fun.getDimension(hds)));
		fun.evaluate(hds, xG, null, G, null);
		System.out.println("Grad: " + G.getVec());
		setHDS(hds);
	}
	
	@Test
	public void testCalculateTriangleAngle() throws Exception {
		double a1 = IsothermicUtility.calculateTriangleAngle(PI/2, 0, PI/4);
		Assert.assertEquals(PI/2, a1, 1E-10);
		double a2 = IsothermicUtility.calculateTriangleAngle(PI/2, 0, -PI/4);
		Assert.assertEquals(PI/2, a2, 1E-10);
		double a3 = IsothermicUtility.calculateTriangleAngle(-PI/2 + EPS, PI/2 - EPS, 0);
		Assert.assertEquals(2*EPS, a3, 1E-10);
		double a4 = IsothermicUtility.calculateTriangleAngle(-PI/2 + EPS, PI/2 - EPS, PI/2);
		Assert.assertEquals(PI - 2*EPS, a4, 1E-10);
		double a5 = IsothermicUtility.calculateTriangleAngle(-PI/8, -3*PI/8, PI/2);
		Assert.assertEquals(PI/4, a5, 1E-10);
		double a6 = IsothermicUtility.calculateTriangleAngle(PI/8, 3*PI/8, PI/2);
		Assert.assertEquals(PI/4, a6, 1E-10);
		double a7 = IsothermicUtility.calculateTriangleAngle(-PI/4, PI/4, PI/2);
		Assert.assertEquals(PI/2, a7, 1E-10);
		double a8 = IsothermicUtility.calculateTriangleAngle(-PI/4, PI/4, 0);
		Assert.assertEquals(PI/2, a8, 1E-10);
	}
	
	@Test
	public void testCalculateCirclePatternRadii() throws Exception {
		Map<CoEdge, Double> betaMap = IsothermicUtility.calculateBetasFromAlphas(hds, alphaMap);
		Map<CoEdge, Double> thetaMap = IsothermicUtility.calculateThetasFromBetas(hds, betaMap);
		Map<CoFace, Double> phiMap = IsothermicUtility.calculatePhisFromBetas(hds, betaMap);
		
		Map<CoFace, Double> rhoMap = IsothermicUtility.calculateCirclePatternRhos(hds, thetaMap, phiMap);
		System.out.println(rhoMap);
		Assert.assertEquals(0.0, rhoMap.get(hds.getFace(0)), 1E-10);
		Assert.assertEquals(0.0, rhoMap.get(hds.getFace(2)), 1E-10);
		Assert.assertEquals(0.6139735886337799, rhoMap.get(hds.getFace(1)), 1E-11);
		
		CPLayoutAdapters layoutAdapters = new CPLayoutAdapters();
		CPLayoutAlgorithm<CoVertex, CoEdge, CoFace>
			layout = new CPLayoutAlgorithm<CoVertex, CoEdge, CoFace>(layoutAdapters, layoutAdapters, rhoMap, thetaMap);
		
		layout.execute(hds);
		
		for (CoVertex vertex : hds.getVertices()) {
			System.out.println(vertex + ": " + Arrays.toString(vertex.T));
		}
	}
	
}