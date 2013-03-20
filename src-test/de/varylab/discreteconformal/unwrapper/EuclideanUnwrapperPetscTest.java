package de.varylab.discreteconformal.unwrapper;

import static java.lang.Math.PI;

import java.util.HashMap;
import java.util.Map;

import org.junit.Assert;
import org.junit.Test;

import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.scene.data.Attribute;
import de.jreality.scene.data.DataList;
import de.jreality.util.Input;
import de.jreality.util.NativePathUtility;
import de.jreality.util.SceneGraphUtility;
import de.jtem.jtao.Tao;

public class EuclideanUnwrapperPetscTest {

	static {
		NativePathUtility.set("native");
		Tao.Initialize();
	}
	
	@Test
	public void testStaticUnwrap() throws Exception {
		ReaderOBJ reader01 = new ReaderOBJ();
		Input in01 = new Input("Obj File 01", getClass().getResourceAsStream("tetraflat.obj"));
		SceneGraphComponent c01 = reader01.read(in01);
		IndexedFaceSet ifs = (IndexedFaceSet)SceneGraphUtility.getFirstGeometry(c01);
		Map<Integer, Double> bdMap = new HashMap<Integer, Double>();
		bdMap.put(1, PI / 3);
		bdMap.put(2, PI / 3 + 0.1);
		bdMap.put(3, PI / 3 - 0.1);
		EuclideanUnwrapperPETSc.unwrap(ifs, bdMap);
		
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
	
}
