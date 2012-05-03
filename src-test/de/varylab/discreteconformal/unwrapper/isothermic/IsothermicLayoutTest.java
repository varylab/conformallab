package de.varylab.discreteconformal.unwrapper.isothermic;

import static java.lang.Math.PI;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class IsothermicLayoutTest {

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
		
		
		AdapterSet aSet = new ConformalAdapterSet();
		IsothermicLayout.doTexLayout(hds, angleMap, aSet);
		
		for (CoVertex tv : hds.getVertices()) {
			System.out.println(Arrays.toString(tv.T));
		}
	}
	
}
