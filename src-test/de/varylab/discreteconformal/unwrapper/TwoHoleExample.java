package de.varylab.discreteconformal.unwrapper;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import de.jreality.geometry.IndexedLineSetFactory;
import de.jreality.geometry.IndexedLineSetUtility;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.plugin.JRViewer;
import de.jreality.reader.Readers;
import de.jreality.scene.Appearance;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.IndexedLineSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.scene.data.Attribute;
import de.jreality.scene.data.StorageModel;
import de.jreality.shader.CommonAttributes;
import de.jreality.util.Input;
import de.jreality.util.SceneGraphUtility;

public class TwoHoleExample {

	public static double[][] verts = {
		{0, 0, 0},
		{3, 0, 0},
		{3, 5, 0},
		{0, 5, 0},
		{1, 1, 0},
		{2, 1, 0},
		{2, 2, 0},
		{1, 2, 0},
		{1, 3, 0},
		{2, 3, 0},
		{2, 4, 0},
		{1, 4, 0}
	};
	public static int[][] indices = {{0,1,2,3, 0},{4,5,6,7, 4}, {8,9,10,11,8}};
	static double  rt = Math.PI/2,
			triangleArea = .005;
	static boolean letter = false, doholes = true;
	public static double[] angles = {rt, rt, rt, rt,  rt, rt,  rt, rt,  rt, rt, rt, rt};
	public static void main(String[] args) {
		// test out a domain with a hole in it.
		IndexedLineSetFactory ilsf = new IndexedLineSetFactory();
		verts = Pn.homogenize(null, verts);
		ilsf.setVertexCount(verts.length);
		ilsf.setVertexCoordinates(verts);
		ilsf.setEdgeCount(indices.length);
		ilsf.setEdgeIndices(indices);
		ilsf.update();
		IndexedLineSet[] AB  = null;
			triangleArea = .05;

		SceneGraphComponent triangulation = null;
		try {
			triangulation = Readers.read(new Input(EuclideanUnwrapProblem.class.getResource("TwoHoleExample.off")));
		} catch (IOException e) {
			e.printStackTrace();
		}
		IndexedFaceSet triang = (IndexedFaceSet) SceneGraphUtility.getFirstGeometry(triangulation);

		JRViewer.display(triang);

		SceneGraphComponent sgc = SceneGraphUtility.createFullSceneGraphComponent("sgc");
		Appearance ap = sgc.getAppearance();
		ap.setAttribute("polygonShader.diffuseColor", Color.white);
		ap.setAttribute(CommonAttributes.FACE_DRAW, false);
		ap.setAttribute(CommonAttributes.EDGE_DRAW, true);
		ap.setAttribute(CommonAttributes.VERTEX_DRAW, false);
//		SimpleTextureFactory stf = new SimpleTextureFactory();
//		stf.setType(TextureType.CHECKERBOARD);
//		stf.setColor(1, Color.black);
//		stf.update();
//		
//		Texture2D tex = TextureUtility.createTexture(ap, "polygonShader", stf.getImageData());
//		Matrix m = MatrixBuilder.euclidean().scale(2).getMatrix();
//		tex.setTextureMatrix(m);
		sgc.setGeometry(triang);
			
		HashMap<Integer, Double> indexToAngle = new HashMap<Integer, Double>();
		for (int i = 0; i<4; ++i) 	{
			indexToAngle.put(i, angles[i]);
		}
		List<Integer> holes = new ArrayList<Integer>();
		holes.add(indices[1][0]);
		holes.add(indices[2][0]);		

		try {
			EuclideanUnwrapperPETSc.unwrapcg(triang, 1, indexToAngle, doholes ? holes : null, 10E-6);
			double[][] tv = triang.getVertexAttributes(Attribute.TEXTURE_COORDINATES).toDoubleArrayArray(null);
			System.err.println("unwrapped verts = "+Rn.toString(tv));
			triang.setVertexAttributes(Attribute.COORDINATES, null);
			triang.setVertexAttributes(Attribute.COORDINATES, triang.getVertexAttributes(Attribute.TEXTURE_COORDINATES));
			} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}

		JRViewer.display(sgc);
		

	}

}
