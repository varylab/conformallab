package de.varylab.discreteconformal.uniformization;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import de.jreality.math.Pn;
import de.jreality.plugin.JRViewer;
import de.jreality.plugin.JRViewer.ContentType;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeLayer;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoTextureDomainPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;
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
		boolean check1 = VisualizationUtility.isInsideFundamentalPolygon(v1, p, 0);
		Assert.assertTrue("v1 is inside", check1);
		
		CoVertex v2 = hds.addNewVertex();
		v2.T = new double[]{-0.25, 0.25, 0, 1};
		boolean check2 = VisualizationUtility.isInsideFundamentalPolygon(v2, p, 0);
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
	
	public static void main(String[] args) {
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
		CoVertex v0 = hds.getVertex(0);
		CoVertex v1 = hds.getVertex(1);
		CoVertex v2 = hds.getVertex(2);
		CoVertex v3 = hds.getVertex(3);
		CoVertex v4 = hds.getVertex(4);
		v0.T = new double[]{0,0,0,1};
		v1.T = new double[]{1,0,0,1};
		v2.T = new double[]{1,1,0,1};
		v3.T = new double[]{0,1,0,1};
		v4.T = new double[]{0.5,0.5,0,1};
		
		HalfedgeInterface hif = new HalfedgeInterface();
		hif.addAdapter(new CoTexturePositionAdapter(), true);
		hif.addAdapter(new CoTextureDomainPositionAdapter(), true);
		JRViewer v = new JRViewer();
		v.registerPlugin(hif);
		v.addBasicUI();
		v.addContentSupport(ContentType.Raw);
		v.addContentUI();
		v.startup();
		
		HalfedgeLayer main = hif.getActiveLayer();
		HalfedgeLayer reglued1 = new HalfedgeLayer(hif);
		HalfedgeLayer reglued2 = new HalfedgeLayer(hif);
		hif.addLayer(reglued1);
		hif.addLayer(reglued2);
		main.set(hds);
		reglued1.setName("Reglued 1");
		reglued2.setName("Reglued 2");
		CoFace reglueFace = e1.getLeftFace();
		Assert.assertTrue(HalfEdgeUtils.isValidSurface(hds, true));
		VisualizationUtility.reglueFace(reglueFace, cutInfo, Pn.EUCLIDEAN);
		reglued1.set(hds);
		Assert.assertTrue(HalfEdgeUtils.isValidSurface(hds, true));
		VisualizationUtility.reglueFace(reglueFace, cutInfo, Pn.EUCLIDEAN);
		reglued2.set(hds);
		Assert.assertTrue(HalfEdgeUtils.isValidSurface(hds, true));
	}
	
}
