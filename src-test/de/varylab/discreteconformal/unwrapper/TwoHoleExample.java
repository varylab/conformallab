package de.varylab.discreteconformal.unwrapper;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.junit.Test;

import de.jreality.math.Rn;
import de.jreality.plugin.JRViewer;
import de.jreality.reader.Readers;
import de.jreality.scene.Appearance;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.scene.data.Attribute;
import de.jreality.shader.CommonAttributes;
import de.jreality.util.Input;
import de.jreality.util.NativePathUtility;
import de.jreality.util.SceneGraphUtility;

public class TwoHoleExample {

	static {
		NativePathUtility.set("native");
	}
	
	static double  rt = Math.PI/2,
			triangleArea = .005;
	static boolean letter = true, doholes = true, originalAngles = false;
	@Test
	public void testLetterB() {

		SceneGraphComponent triangulation = null;
		try {
			triangulation = Readers.read(new Input(EuclideanUnwrapProblem.class.getResource("letterB-TA01.off")));
		} catch (IOException e) {
			e.printStackTrace();
		}
		IndexedFaceSet triang2 = (IndexedFaceSet) SceneGraphUtility.getFirstGeometry(triangulation);
		SceneGraphComponent origSgc = SceneGraphUtility.createFullSceneGraphComponent("sgc");
		origSgc.setGeometry(triang2);

		Appearance ap = origSgc.getAppearance();
		ap.setAttribute("polygonShader.diffuseColor", Color.white);
		ap.setAttribute(CommonAttributes.FACE_DRAW, false);
		ap.setAttribute(CommonAttributes.EDGE_DRAW, true);
		ap.setAttribute(CommonAttributes.VERTEX_DRAW, false);
		HashMap<Integer, Double> indexToAngle = new HashMap<Integer, Double>();
		int[]	inds = new int[]{0,5,45,47};
		for (int i = 0; i<inds.length; ++i) 	{
			indexToAngle.put(inds[i], Math.PI/2);
		}
		List<Integer> holes = new ArrayList<Integer>();
		if (doholes)	{
			holes.add(51);
			holes.add(80);
		}

		try {
			EuclideanUnwrapperPETSc.unwrapcg(triang2, 1, indexToAngle, doholes ? holes : null, 1E-4);
			double[][] tv = triang2.getVertexAttributes(Attribute.TEXTURE_COORDINATES).toDoubleArrayArray(null);
			System.err.println("unwrapped verts = "+Rn.toString(tv));
			triang2.setVertexAttributes(Attribute.COORDINATES, null);
			triang2.setVertexAttributes(Attribute.COORDINATES, triang2.getVertexAttributes(Attribute.TEXTURE_COORDINATES));
			} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}

		JRViewer.display(origSgc);
		

	}

}
