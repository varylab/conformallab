package de.varylab.discreteconformal.uniformization;

import static de.jreality.math.Pn.HYPERBOLIC;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
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

	public CoHDS
		lawson = new CoHDS();
	public AdapterSet
		a = new ConformalAdapterSet();
	public HyperbolicUnwrapper 
		unwrapper = new HyperbolicUnwrapper();
	public FundamentalPolygon 
 		minimalPolygon = null;
	
	@Before
	public void init() throws Exception {
		NativePathUtility.set("native");
		loadSurface();
		unwrapSurface();
		calculateFundamentalPoygon();
	}
	
	
	public void calculateFundamentalPoygon() {
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = unwrapper.getCutInfo();
		FundamentalPolygon cuttedPolygon = FundamentalPolygonUtility.constructFundamentalPolygon(cutInfo, HYPERBOLIC);
		FundamentalVertex root = cuttedPolygon.getMaxValenceVertex();
		minimalPolygon = FundamentalPolygonUtility.minimize(cuttedPolygon, root);
		System.out.println(minimalPolygon);
	}
	
	public void unwrapSurface() throws Exception {
		unwrapper.setGradientTolerance(1E-5);
		unwrapper.setMaxIterations(200);
		unwrapper.unwrap(lawson, 2, a);
	}

	
	public void loadSurface() throws Exception {
		// load lawsons surface
		ReaderOBJ readerOBJ = new ReaderOBJ();
		SceneGraphComponent c = readerOBJ.read(Input.getInput("lawson input", UniformizerTest.class.getResourceAsStream("lawson2498.obj")));
		IndexedFaceSet ifs = (IndexedFaceSet)SceneGraphUtility.getFirstGeometry(c);
		ConverterJR2Heds converter = new ConverterJR2Heds();
		converter.ifs2heds(ifs, lawson, a);
	}
	
	@Test
	public void testCreateSurfaceCurves() throws Exception {
		CoHDS curves = SurfaceCurveUtility.createSurfaceCurves(minimalPolygon, lawson, a, 0, 0.0, true, true);
		Assert.assertEquals("Surface Curve Edge Number", 612, curves.numEdges() / 2);
	}
	
}
