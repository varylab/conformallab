package de.varylab.discreteconformal.heds;


import static de.jreality.shader.CommonAttributes.DIFFUSE_COLOR;
import static de.jreality.shader.CommonAttributes.FACE_DRAW;
import static de.jreality.shader.CommonAttributes.POLYGON_SHADER;
import static de.jreality.shader.CommonAttributes.TUBE_RADIUS;
import static de.jreality.shader.CommonAttributes.TUBE_RADIUS_DEFAULT;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static de.varylab.discreteconformal.util.HomologyUtility.getGeneratorPaths;
import static java.awt.Color.LIGHT_GRAY;

import java.io.IOException;
import java.util.List;
import java.util.Random;
import java.util.Set;

import junit.framework.Assert;

import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.geometry.IndexedFaceSetUtility;
import de.jreality.plugin.JRViewer;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.Appearance;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.Input;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.impl.DoubleArrayAdapter;
import de.jtem.halfedgetools.adapter.type.Color;
import de.jtem.halfedgetools.adapter.type.Radius;
import de.jtem.halfedgetools.jreality.ConverterHeds2JR;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
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
			Input in = new Input("Obj File", HomologyTest.class.getResourceAsStream("brezel2.obj"));
			c =reader.read(in);
			ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
			ConverterJR2Heds converter = new ConverterJR2Heds();
			hds = new CoHDS();
			AdapterSet a = new AdapterSet(new CoPositionAdapter());
			converter.ifs2heds(ifs, hds, a, null);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	@Test
	public void testHomology() throws Exception{
		List<Set<CoEdge>> paths = getGeneratorPaths(hds.getVertex(0), new DefaultWeightAdapter<CoEdge>());
		Assert.assertEquals(4, paths.size());
	}
	
	
	@Color
	private static class EdgeColorAdapter extends DoubleArrayAdapter {

		private double[][]
		    colors = null;
		private double[]
		    defaultColor = {0,0,0};
		private List<Set<CoEdge>> 
			paths = null;
		
		
		public EdgeColorAdapter(List<Set<CoEdge>> paths) {
			super(true, false);
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
		public <
			V extends Vertex<V,E,F>, 
			E extends Edge<V,E,F>, 
			F extends Face<V,E,F>
		> double[] getE(E edge, AdapterSet a) {
			for (Set<CoEdge> path : paths) {
				if (path.contains(edge) || path.contains(edge.getOppositeEdge())) {
					return colors[paths.indexOf(path)];
				}
			}
			return defaultColor;
		}

		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return Edge.class.isAssignableFrom(nodeClass);
		}

		@Override
		public double getPriority() {
			return 0;
		}
		
	}
	
	
	@Radius
	private static class EdgeRadiusAdapter extends AbstractAdapter<Double> {

		private List<Set<CoEdge>> 
			paths = null;
		
		public EdgeRadiusAdapter(List<Set<CoEdge>> paths) {
			super(Double.class, true, false);
			this.paths = paths;
		}
		
		@Override
		public <
			V extends Vertex<V,E,F>, 
			E extends Edge<V,E,F>, 
			F extends Face<V,E,F>
		> Double getE(E edge, AdapterSet a) {
			for (Set<CoEdge> path : paths) {
				if (path.contains(edge) || path.contains(edge.getOppositeEdge())) {
					return 10.0;
				}
			}
			return 0.0;
		}

		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return Edge.class.isAssignableFrom(nodeClass);
		}

		@Override
		public double getPriority() {
			return 0;
		}
		
	}
	
	
	public static void main(String[] args) throws Exception {
		HomologyTest.setUpBeforeClass();
		
		Random rnd = new Random();
		CoVertex root = hds.getVertex(rnd.nextInt(hds.numVertices()));
		
		List<Set<CoEdge>> paths = HomologyUtility.getGeneratorPaths(root, new DefaultWeightAdapter<CoEdge>());
		
		CoPositionAdapter positionAdapter = new CoPositionAdapter();
		EdgeColorAdapter colorAdapter = new EdgeColorAdapter(paths);
		EdgeRadiusAdapter radiusAdapter = new EdgeRadiusAdapter(paths);
		
		ConverterHeds2JR converter = new ConverterHeds2JR();
		AdapterSet a = new AdapterSet();
		a.add(positionAdapter);
		a.add(colorAdapter);
		a.add(radiusAdapter);
		IndexedFaceSet ifs = converter.heds2ifs(hds,  a);
		
		SceneGraphComponent c = new SceneGraphComponent();
		c.setGeometry(ifs);
		IndexedFaceSetUtility.calculateAndSetFaceNormals(ifs);
		Appearance app = new Appearance();
		app.setAttribute(TUBE_RADIUS, TUBE_RADIUS_DEFAULT / 10.0);
		app.setAttribute(VERTEX_DRAW, false);
		app.setAttribute(FACE_DRAW, true);
		app.setAttribute(POLYGON_SHADER + "." + DIFFUSE_COLOR, LIGHT_GRAY);
		c.setAppearance(app);
		JRViewer.display(c);
	}
	
	
	

}
