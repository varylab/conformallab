package de.varylab.discreteconformal.uniformization;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import de.jreality.math.Pn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.math.RnBig;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class VisualizationUtilityTest {

	
	@Test
	public void testIsInsideFundamentalPolygon() throws Exception {
		FundamentalEdge e1 = new FundamentalEdge(0);
		FundamentalEdge e2 = new FundamentalEdge(1);
		FundamentalEdge e3 = new FundamentalEdge(2);
		e1.nextEdge = e2;
		e2.nextEdge = e3;
		e3.nextEdge = e1;
		e1.startPosition = RnBig.toBig(null, new double[]{0,0,0,1});
		e2.startPosition = RnBig.toBig(null, new double[]{1,0,0,1});
		e3.startPosition = RnBig.toBig(null, new double[]{0,1,0,1});
		List<FundamentalEdge> edgeList = new LinkedList<FundamentalEdge>();
		Collections.addAll(edgeList, e1, e2, e3);
		FundamentalPolygon p = new FundamentalPolygon(edgeList);
		
		CoHDS hds = new CoHDS();
		CoVertex v1 = hds.addNewVertex();
		v1.T = new double[]{0.25, 0.25, 0, 1};
		boolean check1 = VisualizationUtility.isInsideFundamentalPolygon(v1, p);
		Assert.assertTrue("v1 is inside", check1);
		
		CoVertex v2 = hds.addNewVertex();
		v2.T = new double[]{-0.25, 0.25, 0, 1};
		boolean check2 = VisualizationUtility.isInsideFundamentalPolygon(v2, p);
		Assert.assertFalse("v2 is outside", check2);
	}
	
	
	@Test
	public void testReglueFace() throws Exception {
		CoHDS hds = new CoHDS();
		CoFace f = HalfEdgeUtils.addNGon(hds, 4);
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
		CoEdge e1 = f.getBoundaryEdge();
		CoEdge e2 = e1.getNextEdge();
		CoEdge e3 = e2.getNextEdge();
		CoEdge e4 = e3.getNextEdge();
		cutInfo.edgeCutMap.put(e1, e3);
		cutInfo.edgeCutMap.put(e1.getOppositeEdge(), e3.getOppositeEdge());
		cutInfo.edgeCutMap.put(e2, e4);
		cutInfo.edgeCutMap.put(e4.getOppositeEdge(), e4.getOppositeEdge());
		cutInfo.cutRoot = e1.getStartVertex();
		TopologyAlgorithms.splitFace(f);
		
		CoFace reglueFace = e1.getLeftFace();
		Assert.assertTrue(HalfEdgeUtils.isValidSurface(hds, true));
		VisualizationUtility.reglueFace(reglueFace, cutInfo, Pn.EUCLIDEAN);
		Assert.assertTrue(HalfEdgeUtils.isValidSurface(hds, true));
		VisualizationUtility.reglueFace(reglueFace, cutInfo, Pn.EUCLIDEAN);
		Assert.assertTrue(HalfEdgeUtils.isValidSurface(hds, true));
	}
	
}
