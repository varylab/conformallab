package de.varylab.discreteconformal.plugin;

import geom3d.Point;

import java.util.ArrayList;
import java.util.List;

import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.CalculatorException;
import de.jtem.halfedgetools.adapter.CalculatorSet;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.jtem.jrworkspace.plugin.PluginInfo;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class EllipticRootGenerator extends AlgorithmPlugin {

	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.Generator;
	}
	
	
	@Override
	public String getAlgorithmName() {
		return "Simple Elliptic Curve";
	}
	
	public < 
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS h, CalculatorSet c, HalfedgeInterface hif) throws CalculatorException {
		CoHDS hds = new CoHDS();
		// first copy
		CoVertex v0a = hds.addNewVertex();
		CoVertex v1a = hds.addNewVertex();
		CoVertex v2a = hds.addNewVertex();
		CoVertex v3a = hds.addNewVertex();
		v0a.getPosition().set(1, 1, 1);
		v1a.getPosition().set(1, -1, -1);
		v2a.getPosition().set(-1, 1, -1);
		v3a.getPosition().set(-1, -1, 1);
		
		List<CoFace> faces = new ArrayList<CoFace>();
		for (int i = 0; i < 8; i++) {
			CoFace f = hds.addNewFace();
			faces.add(f);
			CoEdge e1 = hds.addNewEdge();
			CoEdge e2 = hds.addNewEdge();
			CoEdge e3 = hds.addNewEdge();
			e1.linkNextEdge(e2);
			e2.linkNextEdge(e3);
			e3.linkNextEdge(e1);
			e1.setLeftFace(f);
			e2.setLeftFace(f);
			e3.setLeftFace(f);
			switch (i % 4) {
			case 0:
				e1.setTargetVertex(v0a);
				e2.setTargetVertex(v1a);
				e3.setTargetVertex(v2a);
				break;
			case 1:
				e1.setTargetVertex(v1a);
				e2.setTargetVertex(v3a);
				e3.setTargetVertex(v2a);
				break;
			case 2:
				e1.setTargetVertex(v2a);
				e2.setTargetVertex(v3a);
				e3.setTargetVertex(v0a);
				break;
			case 3:
				e1.setTargetVertex(v3a);
				e2.setTargetVertex(v1a);
				e3.setTargetVertex(v0a);		
				break;
			}
		}
		hds.getEdge(0).linkOppositeEdge(hds.getEdge(6));
		hds.getEdge(1).linkOppositeEdge(hds.getEdge(23));
		hds.getEdge(2).linkOppositeEdge(hds.getEdge(3));
		hds.getEdge(4).linkOppositeEdge(hds.getEdge(10));
		hds.getEdge(5).linkOppositeEdge(hds.getEdge(19));
		hds.getEdge(8).linkOppositeEdge(hds.getEdge(9));
		
		hds.getEdge(12).linkOppositeEdge(hds.getEdge(18));
		hds.getEdge(13).linkOppositeEdge(hds.getEdge(11));
		hds.getEdge(14).linkOppositeEdge(hds.getEdge(15));
		hds.getEdge(16).linkOppositeEdge(hds.getEdge(22));
		hds.getEdge(17).linkOppositeEdge(hds.getEdge(7));
		hds.getEdge(20).linkOppositeEdge(hds.getEdge(21));
		
		for (CoVertex v : hds.getVertices()) {
			double[] pos = v.getPosition().get();
			double l = Rn.euclideanNorm(pos);
			Rn.times(pos, 1 / l, pos);
			v.setTextureCoord(new Point(0, 0, 1));
		}
		hif.set(hds);
	}
	

	@Override
	public PluginInfo getPluginInfo() {
		return new PluginInfo("Root Curve Generator", "Stefan Sechelmann");
	}

}
