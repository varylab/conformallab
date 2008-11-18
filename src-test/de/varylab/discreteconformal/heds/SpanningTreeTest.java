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

import de.jreality.geometry.GeometryUtility;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.Appearance;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.ui.viewerapp.ViewerApp;
import de.jreality.util.Input;
import de.jtem.halfedge.jreality.ConverterHeds2JR;
import de.jtem.halfedge.jreality.ConverterJR2Heds;
import de.jtem.halfedge.jreality.adapter.ColorAdapter2Ifs;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.heds.util.SpanningTreeUtility;

public class SpanningTreeTest {

	private static CHDS 	
		hds = null;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c = null;
		IndexedFaceSet ifs = null;
		try {
			Input in = new Input("Obj File", CLayoutTest.class.getResourceAsStream("brezel.obj"));
			c =reader.read(in);
			ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
			ConverterJR2Heds<CVertex, CEdge, CFace> converter = new ConverterJR2Heds<CVertex, CEdge, CFace>(CVertex.class, CEdge.class, CFace.class);
			hds = new CHDS();
			converter.ifs2heds(ifs, hds, new PositionAdapter());
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public Set<CEdge> makeSpanningTree(Set<CEdge> edges) throws Exception{
		CEdge root = edges.iterator().next();
		Set<CEdge> tree = SpanningTreeUtility.getSpanningTree(edges, root);
		return tree;
	}
	
	
	public Set<CEdge> makeDualSpanningTree(Set<CEdge> edges) throws Exception{
		CEdge root = edges.iterator().next();
		Set<CEdge> tree = SpanningTreeUtility.getDualSpanningTree(edges, root);
		return tree;
	}

	
	private static class EdgeColorAdapter implements ColorAdapter2Ifs<CEdge> {

		private double[]
		    normalColor = {0.0, 0.0, 0.0},
		    markedColor = {9.0, 0.3, 0.2},
		    marked2Color = {0.0, 1.0, 0.2},
		    bothColor = {1.0, 1.0, 0.0};
		private Set<CEdge>
			markedEdges = null,
			markedEdges2 = null;
		
		
		public EdgeColorAdapter(Set<CEdge> marked, Set<CEdge> marked2) {
			this.markedEdges = marked;
			this.markedEdges2 = marked2;
		}
		
		@Override
		public double[] getColor(CEdge edge) {
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
		public AdapterType getAdapterType() {
			return AdapterType.EDGE_ADAPTER;
		}
		
	}
	
	
	public static void main(String[] args) throws Exception {
		SpanningTreeTest.setUpBeforeClass();
		SpanningTreeTest test = new SpanningTreeTest();
		Set<CEdge> edges = new HashSet<CEdge>(hds.getEdges());
		Set<CEdge> tree = test.makeSpanningTree(edges);
		edges.removeAll(tree);
		Set<CEdge> treeDual = test.makeDualSpanningTree(edges);
		PositionAdapter positionAdapter = new PositionAdapter();
		EdgeColorAdapter colorAdapter = new EdgeColorAdapter(tree, treeDual);
		
		ConverterHeds2JR<CVertex, CEdge, CFace> converter = new ConverterHeds2JR<CVertex, CEdge, CFace>();
		IndexedFaceSet ifs = converter.heds2ifs(hds, colorAdapter, positionAdapter);
		
		SceneGraphComponent c = new SceneGraphComponent();
		c.setGeometry(ifs);
		GeometryUtility.calculateAndSetFaceNormals(ifs);
		Appearance app = new Appearance();
		app.setAttribute(TUBE_RADIUS, TUBE_RADIUS_DEFAULT / 10.0);
		app.setAttribute(VERTEX_DRAW, false);
		app.setAttribute(FACE_DRAW, true);
		app.setAttribute(POLYGON_SHADER + "." + DIFFUSE_COLOR, LIGHT_GRAY);
		c.setAppearance(app);
		ViewerApp.display(c);
	}
	
}
