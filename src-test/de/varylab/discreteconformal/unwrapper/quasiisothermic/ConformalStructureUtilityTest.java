package de.varylab.discreteconformal.unwrapper.quasiisothermic;

import static java.lang.Math.PI;
import static java.lang.Math.abs;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class ConformalStructureUtilityTest {

	static {
		NativePathUtility.set("native");
		Tao.Initialize();
	}
	
	private final double
		EPS1 = 1E-2,
		EPS2 = 1E-1;
	public CoHDS
		hds = null;
	public CoVertex
		v0 = null;
	public CoEdge
		e0 = null,
		e1 = null,
		e2 = null,
		e3 = null,
		e4 = null,
		e5 = null,
		e6 = null,
		e7 = null;
	public Map<CoEdge, Double> 
		alphaMap = null;
		
	@Before
	public void init() {
		hds = new CoHDS();
		v0 = hds.addNewVertex();
		CoVertex v1 = hds.addNewVertex();
		CoVertex v2 = hds.addNewVertex();
		CoVertex v3 = hds.addNewVertex();
		CoVertex v4 = hds.addNewVertex();
		
		CoFace f0 = HalfEdgeUtils.constructFaceByVertices(hds, v0, v1, v2).getLeftFace();
		CoFace f1 = HalfEdgeUtils.constructFaceByVertices(hds, v0, v2, v3).getLeftFace();
		CoFace f2 = HalfEdgeUtils.constructFaceByVertices(hds, v0, v3, v4).getLeftFace();
		CoFace f3 = HalfEdgeUtils.constructFaceByVertices(hds, v0, v4, v1).getLeftFace();

		e0 = HalfEdgeUtils.findEdgeBetweenFaces(f1, f0);
		e1 = HalfEdgeUtils.findEdgeBetweenFaces(f2, f1);
		e2 = HalfEdgeUtils.findEdgeBetweenFaces(f3, f2);
		e3 = HalfEdgeUtils.findEdgeBetweenFaces(f0, f3);
		e4 = e0.getNextEdge();
		e5 = e1.getNextEdge();
		e6 = e2.getNextEdge();
		e7 = e3.getNextEdge();
		
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
	}
	
	@Test
	public void testCalculatePseudoConformalStructure() throws Exception {
		Map<CoEdge, Double> lcr = ConformalStructureUtility.calculatePseudoConformalStructure(hds, alphaMap);
		Assert.assertEquals(1.0, lcr.get(e0), 1E-8);
		Assert.assertEquals(1.0, lcr.get(e1), 1E-8);
		Assert.assertEquals(1.0, lcr.get(e2), 1E-8);
		Assert.assertEquals(1.0, lcr.get(e3), 1E-8);
		
		alphaMap.put(e4, -PI/4 + EPS2);
		alphaMap.put(e7,  PI/4 - EPS2);
		lcr = ConformalStructureUtility.calculatePseudoConformalStructure(hds, alphaMap);
		Assert.assertEquals(1.0, lcr.get(e0), 1E-8);
		Assert.assertEquals(1.0948375819248541, lcr.get(e1), 1E-8);
		Assert.assertEquals(1.0, lcr.get(e2), 1E-8);
		Assert.assertEquals(0.9133774876834989, lcr.get(e3), 1E-8);
	}
	
	
	@Test
	public void testCalculateVertexCrossRatioProduct() throws Exception {
		Map<CoEdge, Double> lcr = ConformalStructureUtility.calculatePseudoConformalStructure(hds, alphaMap);
		double p = ConformalStructureUtility.calculateVertexCrossRatioProduct(v0, lcr);
		Assert.assertEquals(1.0, p, 1E-8);
		
		alphaMap.put(e4, -PI/4 + EPS1);
		lcr = ConformalStructureUtility.calculatePseudoConformalStructure(hds, alphaMap);
		p = ConformalStructureUtility.calculateVertexCrossRatioProduct(v0, lcr);
		Assert.assertEquals(1.0202027004321585, p, 1E-8);
	}
	
	
	@Test
	public void testCalculateConformalStructure() throws Exception {
		alphaMap.put(e4, -PI/4 + EPS2);
		Map<CoEdge, Double> lcrPseudoMap = ConformalStructureUtility.calculatePseudoConformalStructure(hds, alphaMap);
		double pPseudo = ConformalStructureUtility.calculateVertexCrossRatioProduct(v0, lcrPseudoMap);
		Assert.assertFalse(abs(pPseudo - 1.0) < 1E-8);
		
		Map<CoEdge, Double> lcr = ConformalStructureUtility.calculateConformalStructure(hds, lcrPseudoMap);
		double p = ConformalStructureUtility.calculateVertexCrossRatioProduct(v0, lcr);
		Assert.assertEquals(1.0, p, 1E-8);
	}
	
	
}
