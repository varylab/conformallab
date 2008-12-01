package de.varylab.discreteconformal.frontend.controller;

import static de.jtem.halfedge.util.HalfEdgeUtils.isValidSurface;

import java.util.LinkedList;
import java.util.List;

import de.jreality.scene.IndexedFaceSet;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.jreality.ConverterHeds2JR;
import de.jtem.halfedge.jreality.ConverterJR2Heds;
import de.jtem.halfedge.jreality.adapter.CoordinateAdapter2Ifs;
import de.jtem.halfedge.jreality.adapter.TextCoordsAdapter2Ifs;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.heds.adapter.TexCoordAdapter;
import de.varylab.discreteconformal.heds.bsp.KdTree;

public class GeometryController {

	private CoHDS
		heds = null;
	private IndexedFaceSet
		ifs = new IndexedFaceSet();
	private KdTree<CoVertex>
		kdTree = null;
	private LinkedList<GeometryChangedListener>
		changeListener = new LinkedList<GeometryChangedListener>();
	private CoordinateAdapter2Ifs<CoVertex>
		coordAdapter = new PositionAdapter();
	private TextCoordsAdapter2Ifs<CoVertex>
		texCoordAdapter = new TexCoordAdapter(true);
	
	public static interface GeometryChangedListener{
		
		public void geometryChanged(CoHDS heds);
		
	}
	
	public GeometryController() {

	}
	
	
	public boolean addChangeListener(GeometryChangedListener l) {
		return changeListener.add(l);
	}
	
	public boolean removeChangeListener(GeometryChangedListener l) {
		return changeListener.remove(l);
	}
	
	private void fireGeometryChanged() {
		for (GeometryChangedListener l : changeListener)
			l.geometryChanged(getCHDS());
	}
	
	public CoHDS getCHDS() {
		return heds;
	}
	
	public IndexedFaceSet getIndexedFaceSet() {
		return ifs;
	}
	
	public KdTree<CoVertex> getKdTree() {
		return kdTree;
	}

	public void setGeometry(CoHDS heds) {
		ConverterHeds2JR<CoVertex, CoEdge, CoFace> converter = new ConverterHeds2JR<CoVertex, CoEdge, CoFace>();
		this.heds = heds;
		if (heds.isTexCoordinatesValid()) {
			ifs = converter.heds2ifs(heds, coordAdapter, texCoordAdapter);
		} else {
			ifs = converter.heds2ifs(heds, coordAdapter);
		}
		generateKdTree();
		fireGeometryChanged();
	}
	
	public void setGeometry(IndexedFaceSet ifs) {
		ConverterJR2Heds<CoVertex, CoEdge, CoFace> converter = new ConverterJR2Heds<CoVertex, CoEdge, CoFace>(CoVertex.class, CoEdge.class, CoFace.class);
		heds = new CoHDS();
		
		converter.ifs2heds(ifs, heds, new PositionAdapter());
		triangulateQuadMesh(heds);
			
		isValidSurface(heds, true);
		this.ifs = ifs;
		generateKdTree();
		fireGeometryChanged();
	}
	
	
	private void generateKdTree() {
		if (heds == null)
			throw new RuntimeException("No Halfedge Datastructure loaded");
		CoVertex[] points = new CoVertex[heds.numVertices()];
		points = heds.getVertices().toArray(points);
		kdTree = new KdTree<CoVertex>(points, 10, false);
	}
	
	
	public static void triangulateQuadMesh(HalfEdgeDataStructure<CoVertex, CoEdge, CoFace> graph) {
    	LinkedList<CoFace> faces = new LinkedList<CoFace>(graph.getFaces());
    	for (CoFace f : faces) {
    		List<CoEdge> b = HalfEdgeUtils.boundaryEdges(f);
    		if (b.size() != 4)
    			continue;
    		CoEdge ea = b.get(0);
    		CoEdge eb = b.get(1);
    		CoEdge ec = b.get(2);
    		CoEdge ed = b.get(3);
    		CoVertex va = ea.getTargetVertex();
    		CoVertex vc = ec.getTargetVertex();
    		CoFace f1 = ea.getLeftFace();
    		CoFace fn = null;
    		if (f1 != null)
    			fn = graph.addNewFace();
    		CoEdge en1 = graph.addNewEdge();
    		CoEdge en2 = graph.addNewEdge();
    		
    		ea.linkNextEdge(en1);
    		ec.linkNextEdge(en2);
    		en1.linkNextEdge(ed);
    		en2.linkNextEdge(eb);
    		
    		ea.setLeftFace(f1);
    		en1.setLeftFace(f1);
    		ed.setLeftFace(f1);
    		
    		eb.setLeftFace(fn);
    		ec.setLeftFace(fn);
    		en2.setLeftFace(fn);
    		
    		en1.setTargetVertex(vc);
    		en2.setTargetVertex(va);
    		
    		en1.linkOppositeEdge(en2);
    	}
    	
    }


	public CoordinateAdapter2Ifs<CoVertex> getCoordAdapter() {
		return coordAdapter;
	}


	public void setCoordAdapter(CoordinateAdapter2Ifs<CoVertex> coordAdapter) {
		this.coordAdapter = coordAdapter;
	}


	public TextCoordsAdapter2Ifs<CoVertex> getTexCoordAdapter() {
		return texCoordAdapter;
	}


	public void setTexCoordAdapter(TextCoordsAdapter2Ifs<CoVertex> texCoordAdapter) {
		this.texCoordAdapter = texCoordAdapter;
	}


}
