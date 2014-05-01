package de.varylab.discreteconformal.plugin.algorithm;

import static de.jreality.geometry.IndexedFaceSetUtility.calculateAndSetVertexNormals;

import java.awt.Color;
import java.util.Set;

import de.jreality.geometry.Primitives;
import de.jreality.math.MatrixBuilder;
import de.jreality.math.Rn;
import de.jreality.scene.Appearance;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.scene.Transformation;
import de.jreality.shader.CommonAttributes;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.jtem.jrworkspace.plugin.PluginInfo;
import de.jtem.numericalMethods.util.Arrays;
import de.varylab.discreteconformal.math.MatrixUtility;
import de.varylab.discreteconformal.plugin.image.ImageHook;

public class AddConeOfRevolutionCommand extends AlgorithmPlugin {

	private IndexedFaceSet
		coneGeometry = Primitives.cone(500);
	private SceneGraphComponent
		cone01 = new SceneGraphComponent(),
		cone02 = new SceneGraphComponent();
	
	public AddConeOfRevolutionCommand() {
		Appearance coneAppearance = new Appearance();
		coneAppearance.setAttribute(CommonAttributes.VERTEX_DRAW, false);
		coneAppearance.setAttribute(CommonAttributes.EDGE_DRAW, false);
		coneAppearance.setAttribute(CommonAttributes.FACE_DRAW, true);
		coneAppearance.setAttribute(CommonAttributes.POLYGON_SHADER + "." + CommonAttributes.DIFFUSE_COLOR, new Color(60, 140, 200));
		coneAppearance.setAttribute(CommonAttributes.POLYGON_SHADER + "." + CommonAttributes.SPECULAR_COLOR, Color.RED);
		coneAppearance.setAttribute(CommonAttributes.POLYGON_SHADER + "." + CommonAttributes.SPECULAR_COEFFICIENT, 1.0);
		coneAppearance.setAttribute(CommonAttributes.POLYGON_SHADER + "." + CommonAttributes.SMOOTH_SHADING, true);
		cone01.setAppearance(coneAppearance);
		cone02.setAppearance(coneAppearance);
		cone01.setGeometry(coneGeometry);
		cone02.setGeometry(coneGeometry);
		cone01.setName("Cone 1");
		cone02.setName("Cone 2");
		MatrixBuilder mb01 = MatrixBuilder.euclidean();
		MatrixBuilder mb02 = MatrixBuilder.euclidean();
		mb01.rotate(Math.PI, 1, 0, 0);
		mb01.translate(0, 0, -1);
		mb02.translate(0, 0, -1);
		mb01.assignTo(cone01);
		mb02.assignTo(cone02);
		calculateAndSetVertexNormals(coneGeometry);
	}
	
	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds, AdapterSet a, HalfedgeInterface hcp) throws Exception {
		Set<V> vSel = hcp.getSelection().getVertices(hds);
		if (vSel.size() == 0) {
			throw new Exception("Select a vertex to define cone size");
		}
		double[] dir = a.getD(Position4d.class, vSel.iterator().next()).clone();
		double[] apex = {0,0,0,1};
		double[] normal = {0,0,1,0};
		Rn.times(dir, 1.2, dir);
		dir[3] = 0.0;
		Transformation coneTrafo = new Transformation("Cone Trafo", getConeMatrix(apex, normal, dir));
		Appearance coneApp = new Appearance("Cone Appearance");
		SceneGraphComponent coneRoot = new SceneGraphComponent("Cone");
		coneRoot.setAppearance(coneApp);
		coneRoot.setTransformation(coneTrafo);
		coneRoot.addChild(cone01);
		coneRoot.addChild(cone02);
		hcp.addTemporaryGeometry(coneRoot);
	}
	
	public double[] getConeMatrix(double[] apex, double[] normal, double[] dir) {
		Rn.normalize(normal, normal);
		double[] ty = Rn.crossProduct(null, normal, dir);
		ty = Arrays.resize(Rn.normalize(ty, ty), 4);
		double[] tx = Rn.crossProduct(null, ty, normal);
		tx = Arrays.resize(Rn.normalize(tx, tx), 4);
		double dot = Rn.innerProduct(tx, dir);
		Rn.times(tx, dot, tx);
		Rn.times(ty, dot, ty);
		double[] e4 = {0,0,0,1};
		double[][] from = {{1,0,0,0}, {0,1,0,0}, {1,0,1,0}, e4};
		double[][] to = {tx, ty, dir, e4};
		double[] T = MatrixUtility.makeMappingMatrix(from, to);
		MatrixBuilder mb = MatrixBuilder.euclidean();
		mb.translate(apex);
		mb.times(T);
		return mb.getMatrix().getArray();
	}

	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.TextureRemeshing;
	}

	@Override
	public String getAlgorithmName() {
		return "Add Cone Of Revolution";
	}

	@Override
	public PluginInfo getPluginInfo() {
		PluginInfo info = super.getPluginInfo();
		info.icon = ImageHook.getIcon("shape_flip_horizontal.png");
		return info;
	}
	
}
