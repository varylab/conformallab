package de.varylab.discreteconformal.frontend.controller;

import static de.jtem.halfedge.util.HalfEdgeUtils.isValidSurface;

import java.util.LinkedList;
import java.util.List;

import de.jreality.scene.IndexedFaceSet;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.jReality.converter.ConverterHeds2JR;
import de.jtem.halfedge.jReality.converter.ConverterJR2Heds;
import de.jtem.halfedge.jReality.interfaces.CoordinateAdapter2Ifs;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.heds.CEdge;
import de.varylab.discreteconformal.heds.CFace;
import de.varylab.discreteconformal.heds.CVertex;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.heds.adapter.TexCoordAdapter;
import de.varylab.discreteconformal.heds.adapter.VertexLabelAdapter;
import de.varylab.discreteconformal.heds.bsp.KdTree;

public class GeometryController {

	private CHDS
		heds = null;
	private IndexedFaceSet
		ifs = new IndexedFaceSet();
	private KdTree<CVertex>
		kdTree = null;
	private LinkedList<GeometryChangedListener>
		changeListener = new LinkedList<GeometryChangedListener>();
	private CoordinateAdapter2Ifs
		coordAdapter = new PositionAdapter();
	
	public static interface GeometryChangedListener{
		
		public void geometryChanged(CHDS heds);
		
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
	
	public CHDS getCHDS() {
		return heds;
	}
	
	public IndexedFaceSet getIndexedFaceSet() {
		return ifs;
	}
	
	public KdTree<CVertex> getKdTree() {
		return kdTree;
	}

	public void setGeometry(CHDS heds) {
		ConverterHeds2JR<CVertex, CEdge, CFace> converter = new ConverterHeds2JR<CVertex, CEdge, CFace>();
		this.heds = heds;
		if (heds.isTexCoordinatesValid()) {
			ifs = converter.heds2ifs(heds, coordAdapter, new VertexLabelAdapter(), new TexCoordAdapter());
		} else {
			ifs = converter.heds2ifs(heds, coordAdapter, new VertexLabelAdapter());
		}
		generateKdTree();
		fireGeometryChanged();
	}
	
	public void setGeometry(IndexedFaceSet ifs) {
		ConverterJR2Heds<CVertex, CEdge, CFace> converter = new ConverterJR2Heds<CVertex, CEdge, CFace>(CVertex.class, CEdge.class, CFace.class);
		heds = new CHDS();
		
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
		CVertex[] points = new CVertex[heds.numVertices()];
		points = heds.getVertices().toArray(points);
		kdTree = new KdTree<CVertex>(points, 10, false);
	}
	
	
	public static void triangulateQuadMesh(HalfEdgeDataStructure<CVertex, CEdge, CFace> graph) {
    	LinkedList<CFace> faces = new LinkedList<CFace>(graph.getFaces());
    	for (CFace f : faces) {
    		List<CEdge> b = HalfEdgeUtils.boundaryEdges(f);
    		if (b.size() != 4)
    			continue;
    		CEdge ea = b.get(0);
    		CEdge eb = b.get(1);
    		CEdge ec = b.get(2);
    		CEdge ed = b.get(3);
    		CVertex va = ea.getTargetVertex();
    		CVertex vc = ec.getTargetVertex();
    		CFace f1 = ea.getLeftFace();
    		CFace fn = null;
    		if (f1 != null)
    			fn = graph.addNewFace();
    		CEdge en1 = graph.addNewEdge();
    		CEdge en2 = graph.addNewEdge();
    		
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


	public CoordinateAdapter2Ifs getCoordAdapter() {
		return coordAdapter;
	}


	public void setCoordAdapter(CoordinateAdapter2Ifs coordAdapter) {
		this.coordAdapter = coordAdapter;
	}


}
