package de.varylab.discreteconformal.unwrapper.quasiisothermic;

import static java.lang.Math.PI;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Assert;
import org.junit.Test;

import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class QuasiIsothermicLayoutTest {

	static {
		NativePathUtility.set("native");
		Tao.Initialize();
	}
	
	
	@Test
	public void testDoLayout() throws Exception {
		CoHDS hds = new CoHDS();
		CoFace f = HalfEdgeUtils.addNGon(hds, 3);
		CoVertex v = TopologyAlgorithms.splitFace(f);

		Map<CoEdge, Double> angleMap = new HashMap<CoEdge, Double>();
		List<CoEdge> eIn = HalfEdgeUtils.incomingEdges(v);
		angleMap.put(eIn.get(0), -3*PI/8);
		angleMap.put(eIn.get(0).getPreviousEdge(), -PI/4);
		angleMap.put(eIn.get(1), 0.0);
		angleMap.put(eIn.get(1).getPreviousEdge(), PI/4);
		angleMap.put(eIn.get(2), 3*PI/8);
		angleMap.put(eIn.get(2).getPreviousEdge(), PI/2);
		
		Set<CoEdge> eSet = new HashSet<CoEdge>(angleMap.keySet());
		for (CoEdge e : eSet) {
			angleMap.put(e.getOppositeEdge(), angleMap.get(e));
		}
		
		Map<CoFace, Double> orientationMap = QuasiisothermicUtility.calculateOrientationFromAlphas(hds, angleMap);
		
		AdapterSet aSet = new ConformalAdapterSet();
		QuasiisothermicLayout.setInitialLength(2.0);
		QuasiisothermicLayout.doTexLayout(hds, angleMap, orientationMap, aSet);
		
		double[] tv0 = {0, 1, 0, 1};
		double[] tv1 = {-1, 0, 0, 1};
		double[] tv2 = {0, -1, 0, 1};
		double[] tv3 = {-0.41421356237309537, 0, 0, 1};
		double[][] tvVec = {tv0, tv1, tv2, tv3};
		for (CoVertex tv : hds.getVertices()) {
			int i = tv.getIndex();
			Assert.assertArrayEquals(tvVec[i], tv.T, 1E-8);
		}
	}
	
}
