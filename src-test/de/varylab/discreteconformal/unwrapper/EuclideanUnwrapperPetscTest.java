package de.varylab.discreteconformal.unwrapper;

import static java.lang.Math.PI;

import java.util.HashMap;
import java.util.Map;

import org.junit.Assert;
import org.junit.Test;

import de.jreality.geometry.IndexedFaceSetUtility;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.plugin.JRViewer;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.scene.data.Attribute;
import de.jreality.scene.data.DataList;
import de.jreality.util.Input;
import de.jreality.util.NativePathUtility;
import de.jreality.util.SceneGraphUtility;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;

public class EuclideanUnwrapperPetscTest {

	private IndexedFaceSet
		ifs = null;
	
	static {
		NativePathUtility.set("native");
		Tao.Initialize();
	}
	
	@Test
	public void testStaticUnwrap() throws Exception {
		ReaderOBJ reader01 = new ReaderOBJ();
		Input in01 = new Input("Obj File 01", getClass().getResourceAsStream("tetraflat.obj"));
		SceneGraphComponent c01 = reader01.read(in01);
		ifs = (IndexedFaceSet)SceneGraphUtility.getFirstGeometry(c01);
		Map<Integer, Double> bdMap = new HashMap<Integer, Double>();
		bdMap.put(1, PI / 3);
		bdMap.put(2, PI / 3 + 0.1);
		bdMap.put(3, PI / 3 - 0.1);
		EuclideanUnwrapperPETSc.unwrap(ifs, 0, bdMap);
		
		DataList texCoords = ifs.getVertexAttributes(Attribute.TEXTURE_COORDINATES);
		double[][] texArr = texCoords.toDoubleArrayArray(new double[ifs.getNumPoints()][4]);
		
		Pn.dehomogenize(texArr, texArr);
		
		double[] t1 = texArr[1];
		double[] t2 = texArr[2];
		double[] t3 = texArr[3];
		double[] v1 = Rn.subtract(null, t1, t2);
		double[] v3 = Rn.subtract(null, t3, t2);
		double angle = Rn.euclideanAngle(v1, v3);
		Assert.assertEquals(PI/3 + 0.1, angle, 1E-10);
	}
	
	
	@Test
	public void testStaticUnwrapBigModel() throws Exception {
		ReaderOBJ reader01 = new ReaderOBJ();
		Input in01 = new Input("Obj File 01", getClass().getResourceAsStream("letterU-01.obj"));
		SceneGraphComponent c01 = reader01.read(in01);
		ifs = (IndexedFaceSet)SceneGraphUtility.getFirstGeometry(c01);
		
		IndexedFaceSetUtility.makeConsistentOrientation(ifs);
		CoHDS hds = new CoHDS();
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		a.add(new CoTexturePositionAdapter());
		ConverterJR2Heds cTo = new ConverterJR2Heds();
		cTo.ifs2heds(ifs, hds, a);
		
		Map<Integer, Double> bdMap = new HashMap<Integer, Double>();
		bdMap.put(476, PI / 2);
		bdMap.put(445, PI / 2);
		bdMap.put(601, PI / 2);
		bdMap.put(574, PI / 2);
		EuclideanUnwrapperPETSc.unwrap(ifs, 300, bdMap);
		
		DataList texCoords = ifs.getVertexAttributes(Attribute.TEXTURE_COORDINATES);
		double[][] texArr = texCoords.toDoubleArrayArray(new double[ifs.getNumPoints()][4]);
		Pn.dehomogenize(texArr, texArr);
		
		double[] t1 = texArr[37];
		double[] t2 = texArr[218];
		double[] t3 = texArr[219];
		double[] v1 = Rn.subtract(null, t1, t2);
		double[] v3 = Rn.subtract(null, t3, t2);
		double angle = Rn.euclideanAngle(v1, v3);
		Assert.assertEquals(PI, angle, 1E-7);
		
		t1 = texArr[397];
		t2 = texArr[445];
		t3 = texArr[427];
		v1 = Rn.subtract(null, t1, t2);
		v3 = Rn.subtract(null, t3, t2);
		angle = Rn.euclideanAngle(v1, v3);
		Assert.assertEquals(PI/2, angle, 1E-7);
		
		Assert.assertEquals(texArr[300][0], 0.0, 1E-7);
		Assert.assertEquals(texArr[300][1], 0.0, 1E-7);
	}
	
	public static void main(String[] args) throws Exception {
		EuclideanUnwrapperPetscTest test = new EuclideanUnwrapperPetscTest();
		test.testStaticUnwrapBigModel();
		JRViewer viewer = new JRViewer();
		viewer.addContentUI();
		viewer.addBasicUI();
		viewer.setContent(test.ifs);
		viewer.startup();
	}
	
}
