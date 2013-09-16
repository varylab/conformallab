package de.varylab.discreteconformal.unwrapper;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import de.jreality.math.MatrixBuilder;
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
	public static int[][] indices  = {{0, 12, 13, 1,2,3, 0},{4,5,6,7, 4}, {8,9,10,11,8}}; //{{0,1,2,3, 0},{4,5,6,7, 4}, {8,9,10,11,8}};
	static double  rt = Math.PI/2,
			triangleArea = .005;
	static boolean letter = true, doholes = true, originalAngles = true;
	public static double[] angles = {rt, rt, rt, rt};
	public static void main(String[] args) {

			triangleArea = .05;

		SceneGraphComponent triangulation = null;
		try {
			triangulation = Readers.read(new Input(EuclideanUnwrapProblem.class.getResource("letterA-03.off")));
		} catch (IOException e) {
			e.printStackTrace();
		}
		IndexedFaceSet triang2 = (IndexedFaceSet) SceneGraphUtility.getFirstGeometry(triangulation);
		SceneGraphComponent origSgc = SceneGraphUtility.createFullSceneGraphComponent("sgc");
		origSgc.setGeometry(triang2);

		try {
			triangulation = Readers.read(new Input(EuclideanUnwrapProblem.class.getResource("letterA-03.off")));
		} catch (IOException e) {
			e.printStackTrace();
		}
		IndexedFaceSet triang = (IndexedFaceSet) SceneGraphUtility.getFirstGeometry(triangulation);
		SceneGraphComponent sgc = SceneGraphUtility.createFullSceneGraphComponent("sgc");
		sgc.setGeometry(triang);

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
		double[] Aangles =  {  1.1899373374161235
				,  1.937769234708022
				,  4.345416072471565
				,  4.363856259874096
				,  1.91932904730549
				,  1.1993786306578298
				,  1.9422140229319633
				,  1.9516553161736694};
	
		HashMap<Integer, Double> indexToAngle = new HashMap<Integer, Double>();
		if (originalAngles) angles = Aangles;
		for (int i = 0; i<angles.length; ++i) 	{
			indexToAngle.put(i, angles[i]);
		}
		List<Integer> holes = new ArrayList<Integer>();
		if (doholes)	{
			holes.add(9);
//			holes.add(indices[2][0]);					
		}

		try {
			EuclideanUnwrapperPETSc.unwrapcg(triang, 1, indexToAngle, doholes ? holes : null, 1E-6);
			double[][] tv = triang.getVertexAttributes(Attribute.TEXTURE_COORDINATES).toDoubleArrayArray(null);
			System.err.println("unwrapped verts = "+Rn.toString(tv));
			triang.setVertexAttributes(Attribute.COORDINATES, null);
			triang.setVertexAttributes(Attribute.COORDINATES, triang.getVertexAttributes(Attribute.TEXTURE_COORDINATES));
			} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		SceneGraphComponent world = SceneGraphUtility.createFullSceneGraphComponent();
		MatrixBuilder.euclidean().translate(2,0,0).assignTo(sgc);
		world.addChildren(sgc, origSgc);

		JRViewer.display(world);
		

	}

}
