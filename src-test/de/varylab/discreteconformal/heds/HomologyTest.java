package de.varylab.discreteconformal.heds;


import static de.jreality.shader.CommonAttributes.DIFFUSE_COLOR;
import static de.jreality.shader.CommonAttributes.FACE_DRAW;
import static de.jreality.shader.CommonAttributes.POLYGON_SHADER;
import static de.jreality.shader.CommonAttributes.TUBE_RADIUS;
import static de.jreality.shader.CommonAttributes.TUBE_RADIUS_DEFAULT;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static java.awt.Color.LIGHT_GRAY;

import java.io.IOException;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.geometry.IndexedFaceSetUtility;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.Appearance;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.ui.viewerapp.ViewerApp;
import de.jreality.util.Input;
import de.jtem.halfedgetools.jreality.ConverterHeds2JR;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.jtem.halfedgetools.jreality.adapter.ColorAdapter2Ifs;
import de.jtem.halfedgetools.jreality.adapter.RelRadiusAdapter2Ifs;
import de.varylab.discreteconformal.adapter.PositionAdapter;
import de.varylab.discreteconformal.util.HomologyUtility;
import de.varylab.discreteconformal.util.Search.DefaultWeightAdapter;

public class HomologyTest {

	private static CoHDS 	
		hds = null;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c = null;
		IndexedFaceSet ifs = null;
		try {
			Input in = new Input("Obj File", EuclideanLayoutTest.class.getResourceAsStream("brezel.obj"));
			c =reader.read(in);
			ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
			ConverterJR2Heds<CoVertex, CoEdge, CoFace> converter = new ConverterJR2Heds<CoVertex, CoEdge, CoFace>(CoVertex.class, CoEdge.class, CoFace.class);
			hds = new CoHDS();
			converter.ifs2heds(ifs, hds, new PositionAdapter());
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	@Test
	public void testHomology() throws Exception{
		System.out.println("HomologyTest.testHomology()");
		System.out.println(hds);
		List<Set<CoEdge>> paths = HomologyUtility.getGeneratorPaths(hds.getVertex(0), new DefaultWeightAdapter<CoEdge>());
		System.out.println("Found " + paths.size() + " generator paths:");
		for (Set<CoEdge> path : paths) {
			System.out.println("Path: length=" + path.size());
		}
	}
	
	
	private static class EdgeColorAdapter implements ColorAdapter2Ifs<CoEdge>, RelRadiusAdapter2Ifs<CoEdge> {

		private double[][]
		    colors = null;
		private double[]
		    defaultColor = {0,0,0};
		private List<Set<CoEdge>> 
			paths = null;
		
		
		public EdgeColorAdapter(List<Set<CoEdge>> paths) {
			this.paths = paths;
			colors = new double[paths.size()][3];
			Random rnd = new Random();
			for (double[] color : colors) {
				color[0] = rnd.nextDouble();
				color[1] = rnd.nextDouble();
				color[2] = rnd.nextDouble();
			}
		}
		
		@Override
		public double[] getColor(CoEdge edge) {
			for (Set<CoEdge> path : paths) {
				if (path.contains(edge) || path.contains(edge.getOppositeEdge())) {
					return colors[paths.indexOf(path)];
				}
			}
			return defaultColor;
		}


		@Override
		public double getReelRadius(CoEdge edge) {
			for (Set<CoEdge> path : paths) {
				if (path.contains(edge) || path.contains(edge.getOppositeEdge())) {
					return 2.0;
				}
			}
			return 0.0;
		}
		
		@Override
		public AdapterType getAdapterType() {
			return AdapterType.EDGE_ADAPTER;
		}
		
	}
	
	
	public static void main(String[] args) throws Exception {
		HomologyTest.setUpBeforeClass();
		
		Random rnd = new Random();
		CoVertex root = hds.getVertex(rnd.nextInt(hds.numVertices()));
		
		List<Set<CoEdge>> paths = HomologyUtility.getGeneratorPaths(root, new DefaultWeightAdapter<CoEdge>());
		
		PositionAdapter positionAdapter = new PositionAdapter();
		EdgeColorAdapter colorAdapter = new EdgeColorAdapter(paths);
		
		ConverterHeds2JR<CoVertex, CoEdge, CoFace> converter = new ConverterHeds2JR<CoVertex, CoEdge, CoFace>();
		IndexedFaceSet ifs = converter.heds2ifs(hds, colorAdapter, positionAdapter);
		
		SceneGraphComponent c = new SceneGraphComponent();
		c.setGeometry(ifs);
		IndexedFaceSetUtility.calculateAndSetFaceNormals(ifs);
		Appearance app = new Appearance();
		app.setAttribute(TUBE_RADIUS, TUBE_RADIUS_DEFAULT / 10.0);
		app.setAttribute(VERTEX_DRAW, false);
		app.setAttribute(FACE_DRAW, true);
		app.setAttribute(POLYGON_SHADER + "." + DIFFUSE_COLOR, LIGHT_GRAY);
		c.setAppearance(app);
		ViewerApp.display(c);
	}
	
	
	

}
