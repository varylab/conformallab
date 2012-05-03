package de.varylab.discreteconformal.unwrapper.isothermic;

import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;
import static java.lang.Math.PI;

import java.util.List;
import java.util.Map;

import junit.framework.Assert;

import org.junit.Test;

import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.jtem.jtao.Tao.Method;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.IsothermicUtility;

public class SinConditionFunctionalTest {

	
	@Test
	public void testSinConditionFunctional() throws Exception {
		NativePathUtility.set("native");
		String[] taoCommand = new String[] {
			"-tao_nm_lamda", "0.01", 
			"-tao_nm_mu", "1.0"
		};
		Tao.Initialize("Sinus Condition Test", taoCommand, false);
		CoHDS hds = new CoHDS();
		CoFace f = HalfEdgeUtils.addNGon(hds, 3);
		CoVertex v = TopologyAlgorithms.splitFace(f);
		Map<Integer, Integer> edgeMap = IsothermicUtility.createUndirectedEdgeMap(hds);
		
		SinConditionFunctional<CoVertex, CoEdge, CoFace, CoHDS> 
			fun = new SinConditionFunctional<CoVertex, CoEdge, CoFace, CoHDS>(hds, edgeMap);
		
		Vec init = new Vec(6);
		List<CoEdge> eIn = HalfEdgeUtils.incomingEdges(v);
		init.setValue(edgeMap.get(eIn.get(0).getIndex()), -3*PI/8, INSERT_VALUES);
		init.setValue(edgeMap.get(eIn.get(0).getPreviousEdge().getIndex()), -PI/4, INSERT_VALUES);
		init.setValue(edgeMap.get(eIn.get(1).getIndex()), 0.2, INSERT_VALUES);
		init.setValue(edgeMap.get(eIn.get(1).getPreviousEdge().getIndex()), PI/4, INSERT_VALUES);
		init.setValue(edgeMap.get(eIn.get(2).getIndex()), 3*PI/8, INSERT_VALUES);
		init.setValue(edgeMap.get(eIn.get(2).getPreviousEdge().getIndex()), PI/2, INSERT_VALUES);
		init.assemble();
		fun.setInitialSolutionVec(init);
		
		Vec startVec = new Vec(init.getSize());
		init.copy(startVec);
		
		Assert.assertTrue("energy is large at the beginning", 1E-3 < fun.evaluateObjective(init));
		
		Tao tao = new Tao(Method.NM);
		tao.setFromOptions();
		tao.setApplication(fun);
		tao.setMaximumIterates(500);
		tao.setTolerances(1E-10, 1E-10, 1E-10, 1E-10);
		tao.solve();
		System.out.println(tao.getSolutionStatus());
		Assert.assertTrue("energy is small after minimizatin", 1E-9 > fun.evaluateObjective(init));
		
		double dif2 = 0.0;
		for (int i = 0; i < startVec.getSize(); i++) {
			double start = startVec.getValue(i);
			double solution = init.getValue(i);
			dif2 += (start-solution)*(start-solution);
		}
		Assert.assertEquals("solution start proximity", 1E-3, dif2, 1E-2);
		
		Tao.Finalize();
	}
	
	
}
