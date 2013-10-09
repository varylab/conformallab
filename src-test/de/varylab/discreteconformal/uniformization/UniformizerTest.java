package de.varylab.discreteconformal.uniformization;

import java.util.Map;

import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.Input;
import de.jreality.util.NativePathUtility;
import de.jreality.util.SceneGraphUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.algorithm.triangulation.Triangulator;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.EuclideanLayout;
import de.varylab.discreteconformal.unwrapper.HyperbolicLayout;

public class UniformizerTest {

	public CoHDS
		wente = null,
		lawson = null;
	public AdapterSet
		a = null;
	
	
	@BeforeClass
	public static void initPetsc() {
		NativePathUtility.set("native");
	}
	
	@Before
	public void loadSurface() throws Exception {
		// load wente torus
		wente = new CoHDS();
		a = new ConformalAdapterSet();
		ReaderOBJ readerOBJ = new ReaderOBJ();
		SceneGraphComponent c = readerOBJ.read(Input.getInput("wente torus input", UniformizerTest.class.getResourceAsStream("wente_torus02.obj")));
		IndexedFaceSet ifs = (IndexedFaceSet)SceneGraphUtility.getFirstGeometry(c);
		ConverterJR2Heds converterJR2Heds = new ConverterJR2Heds();
		converterJR2Heds.ifs2heds(ifs, wente, a);
		Triangulator.triangulateSingleSource(wente);
		
		// load lawsons surface
		lawson = new CoHDS();
		readerOBJ = new ReaderOBJ();
		c = readerOBJ.read(Input.getInput("lawson input", UniformizerTest.class.getResourceAsStream("lawson2498.obj")));
		ifs = (IndexedFaceSet)SceneGraphUtility.getFirstGeometry(c);
		converterJR2Heds.ifs2heds(ifs, lawson, a);
	}
	

	@Test
	public void testGenus1() throws Exception {
		Assert.assertEquals(1240, wente.numVertices(), 1);
 		Map<CoVertex, Double> uMap = Uniformizer.uniformizationFactors(wente, a);
 		Assert.assertEquals(wente.numVertices(), uMap.keySet().size(), 0);
 		for (CoVertex v : wente.getVertices()) {
 			double alpha = EuclideanLayout.calculateAngleSum(v);
 			if (!HalfEdgeUtils.isBoundaryVertex(v)) {
 				Assert.assertEquals(Math.PI * 2, alpha, 1E-5);
 			}
 		}
	}
	
	
	@Test
	public void testGenus2() throws Exception {
		Assert.assertEquals(2498, lawson.numVertices(), 0);
		Map<CoVertex, Double> uMap = Uniformizer.uniformizationFactors(lawson, a);
		Assert.assertEquals(lawson.numVertices(), uMap.keySet().size(), 0);
		for (CoVertex v : lawson.getVertices()) {
 			double alpha = HyperbolicLayout.getAngleSum(v);
 			if (!HalfEdgeUtils.isBoundaryVertex(v)) {
 				Assert.assertEquals(Math.PI * 2, alpha, 1E-5);
 			}
 		}
	}
	
	
}
