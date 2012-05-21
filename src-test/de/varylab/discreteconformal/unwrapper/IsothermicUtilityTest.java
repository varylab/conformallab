package de.varylab.discreteconformal.unwrapper;

import static java.lang.Math.PI;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Assert;
import org.junit.Test;

import cern.colt.Arrays;

import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.IsothermicUtility.CPFunctionalAdapters;
import de.varylab.discreteconformal.unwrapper.IsothermicUtility.CPLayoutAdapters;

public class IsothermicUtilityTest {

	private final double
		EPS = 1E-2;
	
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
		NativePathUtility.set("native");
		Tao.Initialize();
		CoHDS hds = new CoHDS();
		CoFace f = HalfEdgeUtils.addNGon(hds, 3);
		CoVertex v = TopologyAlgorithms.splitFace(f);
		System.out.println("inserted vertex: " + v);

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
		
		CPFunctionalAdapters funAdapters = new CPFunctionalAdapters();
		Vec rho = IsothermicUtility.calculateCirclePatternRadii(hds, angleMap, funAdapters);
		System.out.println(rho);
		Assert.assertEquals(0.0, rho.getValue(0), 1E-10);
		Assert.assertEquals(0.0, rho.getValue(2), 1E-10);
		Assert.assertEquals(0.6139735886337799, rho.getValue(1), 1E-11);
		
		CPLayoutAdapters layoutAdapters = new CPLayoutAdapters(rho);
		CPLayoutAlgorithm<CoVertex, CoEdge, CoFace>
			layout = new CPLayoutAlgorithm<CoVertex, CoEdge, CoFace>(layoutAdapters, layoutAdapters, funAdapters, layoutAdapters, layoutAdapters);
		
		layout.execute(hds);
		
		for (CoVertex vertex : hds.getVertices()) {
			System.out.println(vertex + ": " + Arrays.toString(vertex.T));
		}
	}
	
}
