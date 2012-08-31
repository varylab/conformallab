package de.varylab.discreteconformal.unwrapper.isothermic;

import static java.lang.Math.PI;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import de.jreality.junitutils.Assert;
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

public class IsothermicLayoutTest {

	static {
		NativePathUtility.set("native");
		String[] args = {
				"-help",
				"-snes_view",
				"-snes_type", "ls",
				"-snes_test_display",
				"-ksp_converged_reason",
				"-pc_factor_shift_nonzero", "1.0e-10"
		};
		Tao.Initialize("Sinus Condition Test", args, false);
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
		
		Map<CoFace, Double> orientationMap = IsothermicUtility.calculateOrientationFromAlphas(hds, angleMap);
		
		AdapterSet aSet = new ConformalAdapterSet();
		IsothermicLayout.setInitialLength(2.0);
		IsothermicLayout.doTexLayout(hds, angleMap, orientationMap, aSet);
		
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
