package de.varylab.discreteconformal.heds;


import static de.jreality.shader.CommonAttributes.DIFFUSE_COLOR;
import static de.jreality.shader.CommonAttributes.FACE_DRAW;
import static de.jreality.shader.CommonAttributes.POLYGON_SHADER;
import static de.jreality.shader.CommonAttributes.TUBE_RADIUS;
import static de.jreality.shader.CommonAttributes.TUBE_RADIUS_DEFAULT;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;
import static java.awt.Color.LIGHT_GRAY;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import org.junit.BeforeClass;

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
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.impl.DoubleArrayAdapter;
import de.jtem.halfedgetools.adapter.type.Color;
import de.jtem.halfedgetools.jreality.ConverterHeds2JR;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.util.SpanningTreeUtility;

public class SpanningTreeTest {

	private static CoHDS 	
		hds = null;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c = null;
		IndexedFaceSet ifs = null;
		try {
			Input in = new Input("Obj File", SpanningTreeTest.class.getResourceAsStream("brezel.obj"));
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
	
	
	public Set<CoEdge> makeSpanningTree(Set<CoEdge> edges) throws Exception{
		CoEdge root = edges.iterator().next();
		Set<CoEdge> tree = SpanningTreeUtility.getSpanningTree(edges, root);
		return tree;
	}
	
	
	public Set<CoEdge> makeDualSpanningTree(Set<CoEdge> edges) throws Exception{
		CoEdge root = edges.iterator().next();
		Set<CoEdge> tree = SpanningTreeUtility.getDualSpanningTree(edges, root);
		return tree;
	}

	
	@Color
	private static class EdgeColorAdapter extends DoubleArrayAdapter {

		private double[]
		    normalColor = {0.0, 0.0, 0.0},
		    markedColor = {9.0, 0.3, 0.2},
		    marked2Color = {0.0, 1.0, 0.2},
		    bothColor = {1.0, 1.0, 0.0};
		private Set<CoEdge>
			markedEdges = null,
			markedEdges2 = null;
		
		
		public EdgeColorAdapter(Set<CoEdge> marked, Set<CoEdge> marked2) {
			super(true, false);
			this.markedEdges = marked;
			this.markedEdges2 = marked2;
		}
		
		@Override
		public <
			V extends Vertex<V, E, F>,
			E extends Edge<V, E, F>,
			F extends Face<V, E, F>
		> double[] getE(E edge, AdapterSet a) {
			if (markedEdges.contains(edge) && markedEdges2.contains(edge)) {
				return bothColor;
			} else if (markedEdges.contains(edge)) {
				return markedColor;
			} else if (markedEdges2.contains(edge)){
				return marked2Color;
			} else {
				return normalColor;
			}
		}

		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return CoEdge.class.isAssignableFrom(nodeClass);
		}

		@Override
		public double getPriority() {
			return 0;
		}

	}
	
	
	public static void main(String[] args) throws Exception {
		SpanningTreeTest.setUpBeforeClass();
		SpanningTreeTest test = new SpanningTreeTest();
		Set<CoEdge> edges = new HashSet<CoEdge>(hds.getEdges());
		Set<CoEdge> tree = test.makeSpanningTree(edges);
		edges.removeAll(tree);
		Set<CoEdge> treeDual = test.makeDualSpanningTree(edges);
		CoPositionAdapter positionAdapter = new CoPositionAdapter();
		EdgeColorAdapter colorAdapter = new EdgeColorAdapter(tree, treeDual);
		
		ConverterHeds2JR converter = new ConverterHeds2JR();
		AdapterSet a = new AdapterSet();
		a.add(positionAdapter);
		a.add(colorAdapter);
		IndexedFaceSet ifs = converter.heds2ifs(hds, a);
		
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
