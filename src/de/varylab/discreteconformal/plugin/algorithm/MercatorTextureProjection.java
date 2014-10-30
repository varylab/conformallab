package de.varylab.discreteconformal.plugin.algorithm;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition2d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class MercatorTextureProjection extends StereographicTextureProjection {

	public MercatorTextureProjection(HalfedgeInterface master) {
		super(master);
	}

	@Override
	public String getAlgorithmName() {
		return "Mercator Projection";
	}

	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds2, AdapterSet a, HalfedgeInterface hi) {
		super.execute(hds2, a, hi);
		CoHDS hds = hi.get(new CoHDS());

//		// create branch cut	
//		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
//		double snapTolerance = 1E-20;
//		int signature = Pn.EUCLIDEAN;
//		Set<CoVertex> newVertices = new HashSet<>();
//		double[][] segment = {{-1E10, 0, 1},{0, 0, 1}};
//		createIntersectionVertices(segment, true, hds, hds, cutInfo, snapTolerance, signature, newVertices);
//		triangulateByCuttingCorners(hds, a);
//		Selection pathSelection = new Selection();
//		pathSelection.addAll(newVertices);
//		LinkedList<CoEdge> newCut = DiscreteConformalPlugin.selectCutPath(hds, newVertices, pathSelection);
//		CuttingInfo<CoVertex, CoEdge, CoFace> newCutInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
//		CuttingUtility.cutAlongPath(newCut, newCutInfo);
//		
//		for (CoVertex v : newVertices) {
//			boolean isUpper = false;
//			double[] p = a.getD(TexturePosition2d.class, v);
//			for (CoVertex vv : HalfEdgeUtilsExtra.getVertexStar(v)) {
//				if (HalfEdgeUtils.isBoundaryVertex(vv)) continue;
//				double[] pp = a.getD(TexturePosition2d.class, vv);
//				if (pp[1] > 1E-8) {
//					isUpper = true;
//					break;
//				}
//			}
//			double offset = 1E-8;
//			if (isUpper) {
//				a.set(TexturePosition.class, v, new double[]{p[0], offset});
//				if (newCutInfo.vertexCopyMap.containsKey(v)) {
//					CoVertex vv = newCutInfo.vertexCopyMap.get(v);
//					a.set(TexturePosition.class, vv, new double[]{p[0], -offset});
//				}
//			} else {
//				a.set(TexturePosition.class, v, new double[]{p[0], -offset});
//				if (newCutInfo.vertexCopyMap.containsKey(v)) {
//					CoVertex vv = newCutInfo.vertexCopyMap.get(v);
//					a.set(TexturePosition.class, vv, new double[]{p[0], offset});
//				}
//			}
//		}
		
		// complex logarithm
		for (CoVertex v : hds.getVertices()) {
			double[] pos = a.getD(TexturePosition2d.class, v);
			Complex c = new Complex(pos[0], pos[1]);
			c = c.log();
			a.set(TexturePosition.class, v, new double[]{c.re, c.im});
		}
		hi.update();
		masterInterface.update();
//		masterInterface.setSelection(new Selection(newCutInfo.edgeCutMap.keySet()));
	}

}
