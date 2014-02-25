package de.varylab.discreteconformal.plugin.schottky;

import java.io.InputStream;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.junit.Assert;
import org.junit.Test;

import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.mfc.field.Complex;
import de.varylab.conformallab.data.DataIO;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.SchottkyData;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.NodeIndexComparator;

public class SchottkyUtilityTest {

	static {
		NativePathUtility.set("native");
	}
	
	@Test
	public void testGenerateSurface() throws Exception {
		InputStream dataIn = getClass().getResourceAsStream("genus2.xml");
		SchottkyData sd = DataIO.readConformalData(SchottkyData.class, dataIn);
		List<SchottkyGenerator> data = DataUtility.toSchottkyGeneratorsList(sd);
		CoHDS hds = new CoHDS();
		Map<CoEdge, Double> lMap = new HashMap<CoEdge, Double>();
		List<Set<CoEdge>> cycles = new LinkedList<Set<CoEdge>>();
		Map<CoVertex, double[]> mapCycleMap = new HashMap<CoVertex, double[]>();
		Complex rootPos = new Complex(0, 0);
		SchottkyUtility.generateSurface(hds, data, rootPos, lMap, cycles, mapCycleMap, 0, 10, 10, 10);
		Assert.assertEquals(2, HalfEdgeUtils.getGenus(hds));
	}
	
	
	@Test
	public void testUnwrapSchottkySurface() throws Exception {
		InputStream dataIn = getClass().getResourceAsStream("genus2.xml");
		SchottkyData sd = DataIO.readConformalData(SchottkyData.class, dataIn);
		List<SchottkyGenerator> data = DataUtility.toSchottkyGeneratorsList(sd);
		CoHDS hds = new CoHDS();
		Map<CoEdge, Double> lMap = new HashMap<CoEdge, Double>();
		List<Set<CoEdge>> cycles = new LinkedList<Set<CoEdge>>();
		Map<CoVertex, double[]> mapCycleMap = new HashMap<CoVertex, double[]>();
		Complex rootPos = new Complex(0, 0);
		CoVertex root = SchottkyUtility.generateSurface(hds, data, rootPos, lMap, cycles, mapCycleMap, 0, 4, 0, 0);
		SchottkyLengthAdapter schottkyMetric = new SchottkyLengthAdapter(lMap);
		AdapterSet aSet = new ConformalAdapterSet();
		aSet.add(schottkyMetric);
		
		SchottkyUtility.unwrapSchottkySurface(hds, cycles, mapCycleMap, root, aSet, false, false);
		
		TreeSet<CoVertex> boundary = new TreeSet<CoVertex>(new NodeIndexComparator<CoVertex>());
		boundary.addAll(HalfEdgeUtils.boundaryVertices(hds));
		System.out.println(boundary);
	}
	
	
}
