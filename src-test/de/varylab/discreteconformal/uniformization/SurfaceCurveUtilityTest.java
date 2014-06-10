package de.varylab.discreteconformal.uniformization;

import static de.jreality.math.Pn.HYPERBOLIC;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import de.jreality.math.Pn;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.IndexedLineSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.scene.data.Attribute;
import de.jreality.util.Input;
import de.jreality.util.NativePathUtility;
import de.jreality.util.SceneGraphUtility;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.HyperbolicUnwrapper;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class SurfaceCurveUtilityTest {

	static {
		NativePathUtility.set("native");
	}
	
	public FundamentalPolygon calculateFundamentalPoygon(CoHDS lawson) throws Exception {
		AdapterSet a = new ConformalAdapterSet();
		HyperbolicUnwrapper unwrapper = new HyperbolicUnwrapper();
		unwrapper.setGradientTolerance(1E-5);
		unwrapper.setMaxIterations(200);
		unwrapper.unwrap(lawson, 2, a);
		// calculate fundamental polygon
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = unwrapper.getCutInfo();
		FundamentalPolygon cuttedPolygon = FundamentalPolygonUtility.constructFundamentalPolygon(cutInfo, HYPERBOLIC);
		FundamentalVertex root = cuttedPolygon.getMaxValenceVertex();
		return FundamentalPolygonUtility.minimize(cuttedPolygon, root);
	}
	
	public CoHDS loadSurface() throws Exception {
		AdapterSet a = new ConformalAdapterSet();
		ReaderOBJ readerOBJ = new ReaderOBJ();
		SceneGraphComponent c = readerOBJ.read(Input.getInput("lawson input", UniformizerTest.class.getResourceAsStream("lawson2498.obj")));
		IndexedFaceSet ifs = (IndexedFaceSet)SceneGraphUtility.getFirstGeometry(c);
		ConverterJR2Heds converter = new ConverterJR2Heds();
		CoHDS lawson = new CoHDS();
		converter.ifs2heds(ifs, lawson, a);
		return lawson;
	}
	
	@Test
	public void testCreateSurfaceCurves() throws Exception {
		AdapterSet a = new ConformalAdapterSet();
		CoHDS lawson = loadSurface();
		FundamentalPolygon minimalPolygon = calculateFundamentalPoygon(lawson);
		IndexedLineSet curves = SurfaceCurveUtility.createSurfaceCurves(minimalPolygon, lawson, a, 0, 0.0, true, true, Pn.HYPERBOLIC);
		Assert.assertEquals("Surface Curve Edge Number", 607, curves.getEdgeAttributes(Attribute.INDICES).toIntArrayArray().size());
	}
	
	
	@Test
	public void testIsBetween() {
		double[] x1 = {1,1,1,1};
		double[] x2 = {1 + 1E-5,1,1,1};
		double[][] s = {{0,0,0,1},{2,2,2,1}};
		Assert.assertTrue(SurfaceCurveUtility.isOnSegment(x1, s));
		Assert.assertFalse(SurfaceCurveUtility.isOnSegment(x2, s));
	}
	
	@Test
	public void testGetPointOnSegment_SegmentEdge() {
		double[] edgePoint = {0.352392439203295, 0.9123804930829212, 1.0};
		double[][] edgeSegment = {
			{0.34745306897719913, 0.912568467121888, 1.0}, 
			{0.35239243920431296, 0.912380493082885, 1.0}
		};
		double[][] targetSegment = {
			{-0.2896352574166635, 0.03146361746587523, 0.10643898885661185, 0.44373020886051084}, 
			{-0.2666822290964323, 0.019034256494171405, 0.10525293907970201, 0.44373020886051084}
		}; 
		double[] result = SurfaceCurveUtility.getPointOnCorrespondingSegment(edgePoint, edgeSegment, targetSegment, Pn.HYPERBOLIC);
		for (int i = 0; i < result.length; i++) {
			Assert.assertEquals(targetSegment[1][i], result[i], 1E-10);
		}
	}
	
	@Test@Ignore
	public void testPnDistanceBetweenHyperbolic() {
		double[] p1 = {0.35239243920329500, 0.9123804930829212, 1.0};
		double[] p2 = {0.35239243920431296, 0.9123804930828850, 1.0};
		double dist = Pn.distanceBetween(p1, p2, HYPERBOLIC);
		Assert.assertTrue(!Double.isNaN(dist));
	}
	
}
