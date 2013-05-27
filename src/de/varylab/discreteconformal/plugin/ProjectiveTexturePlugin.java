package de.varylab.discreteconformal.plugin;

import static de.varylab.discreteconformal.unwrapper.quasiisothermic.QuasiisothermicUtility.calculateQuasiconformalFactors;
import static java.lang.Math.exp;

import java.util.HashMap;
import java.util.Map;

import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition2d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.jtem.jrworkspace.plugin.Controller;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class ProjectiveTexturePlugin extends AlgorithmPlugin {

	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.Geometry;
	}
	
	@Override
	public String getAlgorithmName() {
		return "Optimize Texture Interpolation";
	}
	
	@Override
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds, AdapterSet a, HalfedgeInterface hi) {
		CoHDS coHDS = hi.get(new CoHDS());
		calculateProjectiveTexCoords(coHDS, a);
		hi.set(coHDS);
	}
	
	public void calculateProjectiveTexCoords(CoHDS hds, AdapterSet a) {
		Map<CoEdge, Double> edgeLengthMap = new HashMap<CoEdge, Double>();
		Map<CoEdge, Double> texLengthMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : hds.getPositiveEdges()) {
			Double l = a.get(Length.class, e, Double.class);
			double[] t1 = a.getD(TexturePosition2d.class, e.getStartVertex());
			double[] t2 = a.getD(TexturePosition2d.class, e.getTargetVertex());
			double l2 = Rn.euclideanDistance(t1, t2);
			edgeLengthMap.put(e, l);
			edgeLengthMap.put(e.getOppositeEdge(), l);
			texLengthMap.put(e, l2);
			texLengthMap.put(e.getOppositeEdge(), l2);
		}
		Map<CoVertex, Double> uMap = calculateQuasiconformalFactors(hds, edgeLengthMap, texLengthMap);
		for (CoVertex v : hds.getVertices()) {
			Pn.dehomogenize(v.T, v.T);
			double ui = uMap.get(v);
			double e = exp(-ui/2);
			Rn.times(v.T, e, v.T);
		}
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
	}
	
}
