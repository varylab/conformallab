package de.varylab.discreteconformal.plugin.algorithm;

import static de.jreality.math.Pn.EUCLIDEAN;

import java.util.Iterator;
import java.util.Set;

import de.jreality.math.FactoredMatrix;
import de.jreality.math.P2;
import de.jreality.math.P3;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.Adapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeLayer;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.jtem.jrworkspace.plugin.PluginInfo;
import de.varylab.discreteconformal.adapter.ConeMapAdapter;
import de.varylab.discreteconformal.adapter.CylinderMapAdapter;
import de.varylab.discreteconformal.plugin.image.ImageHook;

public class MapToConeCommand extends AlgorithmPlugin {
	
	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS hds, AdapterSet a, HalfedgeInterface hi) throws Exception {
		 Set<E> eSel = hi.getSelection().getEdges(hds);
		 if (eSel.size() != 4) {
			 throw new Exception("Select 2 boundary edges to define cone identification");
		 }
		 Iterator<E> it = eSel.iterator();
		 E e1 = it.next();
		 E e2 = it.next();
		 if (e2 == e1.getOppositeEdge()) {
			 e2 = it.next();
		 }
		 if (e1.getLeftFace() == null) {
			 e1 = e1.getOppositeEdge();
		 }
		 if (e2.getLeftFace() == null) {
			 e2 = e2.getOppositeEdge();
		 }
		 Adapter<double[]> cone = mapToCone(e1, e2, a);
		 HalfedgeLayer l = hi.createLayer("Cone");
		 l.addAdapter(cone, true);
		 l.set(hds);
	}
	
	private <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	>  Adapter<double[]> mapToCone(E ce, E opce, AdapterSet a) throws Exception {
		AdapterSet result = new AdapterSet();
		V v1 = ce.getStartVertex(); 
		V v2 = ce.getTargetVertex();
		V w1 = opce.getTargetVertex();
		V w2 = opce.getStartVertex();
		System.out.println(v1 + "->" + w1);
		System.out.println(v2 + "->" + w2);
		double[]
			p1 = P2.projectP3ToP2(null, a.getD(TexturePosition4d.class, v1)),
			p2 = P2.projectP3ToP2(null, a.getD(TexturePosition4d.class, v2)),
			q1 = P2.projectP3ToP2(null, a.getD(TexturePosition4d.class, w1)),
			q2 = P2.projectP3ToP2(null, a.getD(TexturePosition4d.class, w2));
		
		double[] v1t = a.getD(TexturePosition4d.class, v1);
		double[] v2t = a.getD(TexturePosition4d.class, v2);
		double[] w1t = a.getD(TexturePosition4d.class, w1);
		double[] w2t = a.getD(TexturePosition4d.class, w2);
		
		FactoredMatrix M = new FactoredMatrix(P3.makeDirectIsometryFromFrames(null, v1t, Rn.subtract(null, v2t, v1t), new double[]{0,0,1,0}, w1t, Rn.subtract(null, w2t, w1t), new double[]{0,0,1,0}, EUCLIDEAN));
		if(M.getRotationAngle() == 0) {
			int direction = 0;
			double[] t = M.getTranslation();
			if(Math.abs(t[0]) < 1E-8) {
				direction = 1;
			} else if(Math.abs(t[1]) < 1E-8) {
				direction = 0;
			} else {
				throw new Exception("Texture does not seem to be periodic in horizontal or vertical direction.");
			}
			double dist = Math.abs(t[direction]);
			CylinderMapAdapter cma = a.query(CylinderMapAdapter.class);
			if(cma == null) {
				cma = new CylinderMapAdapter(v1, direction, dist);
				result.add(cma);
			} else {
				cma.initialize(v1,direction,dist);
			}
			return cma;
		} else {
			double[] origin = Pn.dehomogenize(null, P2.pointFromLines(null, 
				perpendicularBisector(p1, q1),
				perpendicularBisector(p2, q2)
			));
			double period = getRotationAngle(p1,q1,origin);
			ConeMapAdapter coneAdapter = a.query(ConeMapAdapter.class);	
			if(coneAdapter == null) {
				coneAdapter = new ConeMapAdapter(origin, period);
				result.add(coneAdapter);
			} else {
				coneAdapter.initialize(origin, period);
			}
			return coneAdapter;
		}
	}

	private double getRotationAngle(double[] p, double[] q, double[] o) {
		p = Pn.dehomogenize(null, p);
		q = Pn.dehomogenize(null, q);
		o = Pn.dehomogenize(null, o);
		double[] v1 = Rn.subtract(null, p, o);
		double[] v2 = Rn.subtract(null, q, o);
		double alphaP = Math.atan2(v1[1], v1[0]);
		double alphaQ = Math.atan2(v2[1], v2[0]);
		return alphaQ - alphaP;
	}

	private double[] perpendicularBisector(double[] p,	double[] q) {
		p = Pn.dehomogenize(null, p);
		q = Pn.dehomogenize(null, q);
		double[] midpoint = Rn.linearCombination(null, 0.5, p, 0.5, q);
		double[] v = Rn.subtract(null, p, q);
		return P2.lineFromPoints(null, midpoint, new double[]{-v[1],v[0],0});
	}

	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.TextureRemeshing;
	}

	@Override
	public String getAlgorithmName() {
		return "Map To Cone";
	}

	@Override
	public PluginInfo getPluginInfo() {
		PluginInfo info = super.getPluginInfo();
		info.icon = ImageHook.getIcon("shape_flip_horizontal.png");
		return info;
	}

}
