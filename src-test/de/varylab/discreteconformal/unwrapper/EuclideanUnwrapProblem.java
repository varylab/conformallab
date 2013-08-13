package de.varylab.discreteconformal.unwrapper;

import java.io.IOException;
import java.util.HashMap;

import de.jreality.geometry.PointSetFactory;
import de.jreality.geometry.PointSetUtility;
import de.jreality.geometry.Primitives;
import de.jreality.math.Matrix;
import de.jreality.math.MatrixBuilder;
import de.jreality.plugin.JRViewer;
import de.jreality.reader.Readers;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.scene.data.Attribute;
import de.jreality.shader.CommonAttributes;
import de.jreality.util.Input;
import de.jreality.util.SceneGraphUtility;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapperPETSc;

public class EuclideanUnwrapProblem {

	static HashMap<Integer, Double> indexToAngle = new HashMap<Integer, Double>();
	static int[] corners = {24, 39, 57, 7, 34};

	public static void main(String[] args) {
		SceneGraphComponent triangulation = null;
		try {
			triangulation = Readers.read(new Input(EuclideanUnwrapProblem.class.getResource("diskUnwrapExample.off")));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		IndexedFaceSet ifs = (IndexedFaceSet) SceneGraphUtility.getFirstGeometry(triangulation);
		int n = ifs.getNumPoints();
		for (int i = 0; i<n ; ++i)	{
			indexToAngle.put(i,Math.PI);
//			System.err.println("angle = "+total);
		}
		for (int i = 0; i<4; ++i)	{
			indexToAngle.put(corners[i],Math.PI/2);
		}
		EuclideanUnwrapperPETSc.unwrap(ifs, corners[4], indexToAngle, null);

		double[][] verts = ifs.getVertexAttributes(Attribute.COORDINATES).toDoubleArrayArray(null);
		double[][] cornerpts = {verts[corners[0]], verts[corners[1]], verts[corners[2]], verts[corners[3]]};
		PointSetFactory psf = new PointSetFactory();
		psf.setVertexCount(4);
		psf.setVertexCoordinates(cornerpts);
		psf.update();
		SceneGraphComponent cornerSGC = SceneGraphUtility.createFullSceneGraphComponent("corners");
		triangulation.addChild(cornerSGC);
		cornerSGC.setGeometry(psf.getGeometry());
		cornerSGC.getAppearance().setAttribute(CommonAttributes.VERTEX_DRAW, true);
		cornerSGC.getAppearance().setAttribute(CommonAttributes.POINT_RADIUS, .2);
		triangulation.getAppearance().setAttribute(CommonAttributes.FACE_DRAW, false);
		JRViewer.display(triangulation);
	}
}
