package de.varylab.discreteconformal.adapter;

import static de.jreality.shader.CommonAttributes.EDGE_DRAW;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;

import java.util.Iterator;
import java.util.LinkedList;

import de.jreality.geometry.IndexedFaceSetUtility;
import de.jreality.math.FactoredMatrix;
import de.jreality.math.P2;
import de.jreality.math.P3;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.plugin.JRViewer;
import de.jreality.scene.Appearance;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.shader.CommonAttributes;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition2d;

public class ConeMapAdapter extends AbstractAdapter<double[]> {

	private double[] origin = null;
	private double period = 1.0;
	
	public ConeMapAdapter(double[] origin, double anglePeriod) {
		super(double[].class, true, false);
		initialize(origin, anglePeriod);
	}
	
	@Override
	public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
		return Vertex.class.isAssignableFrom(nodeClass);
	}

	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double[] getV(V v, AdapterSet a) {
		double[] vT = a.getD(TexturePosition2d.class, v);
		return cone(origin, period, vT);
	}

	private static double[] cone(double[] o, double anglePeriod, double[] vT) {
		double[] pt = Rn.subtract(null, vT, o);
		double alpha = 2*Math.PI*((Math.atan2(pt[1], pt[0])+2*Math.PI) % anglePeriod)/anglePeriod;
		double r = Rn.euclideanNorm(pt);
		double x = 0.5*anglePeriod/Math.PI;
		double[] g = new double[]{ x, Math.sqrt(1.0-x*x)};
		return new double[]{ r*Math.cos(alpha)*g[0], r*Math.sin(alpha)*g[0], r*g[1] };
	}

	public void initialize(double[] origin, double anglePeriod) {
		if(origin.length == 3) {
			origin[0]/=origin[2];
			origin[1]/=origin[2];
		}
		this.origin = new double[]{origin[0], origin[1]};
		period = Math.abs(anglePeriod);
	}
	
	public double[][] getArc(Edge<?,?,?> e, AdapterSet a, int n) {
		
		double[] start = a.getD(TexturePosition2d.class, e.getStartVertex());
		double[] target = a.getD(TexturePosition2d.class, e.getTargetVertex());
		
		double[][] curve = conicalArc(origin, period, n, start, target);
		return curve;
	}

	private static double[][] conicalArc(double[] o, double period, int n, double[] start, double[] target) {
		double[] dif = Rn.subtract(null, target, start);
		double[] step = Rn.times(null, 1.0/n, dif);
		
		LinkedList<double[]> domainVerts = new LinkedList<double[]>();
		
		for(int i = 0; i <= n; ++i) {
			double[] pt = Rn.add(null, start, Rn.times(null, i, step));
			domainVerts.add(pt);
		}
		double[][] curve = new double[domainVerts.size()][3];
		int i = 0;
		for (Iterator<double[]> iterator = domainVerts.iterator(); iterator.hasNext();) {
			curve[i++] = cone(o, period, iterator.next());
		}
		return curve;
	}
	
	public static void main(String[] args) {
		SceneGraphComponent root = new SceneGraphComponent();
		Appearance app = new Appearance();
		app.setAttribute(VERTEX_DRAW, true);
		app.setAttribute(EDGE_DRAW, false);
		app.setAttribute(CommonAttributes.FACE_DRAW, false);
		root.setAppearance(app);
		double[][] arc = new double[][]{ {1.5,1.0,0.0,1.0},{2.0, 1.0,0.0,1.0}, {1.0,2.0,0.0,1.0}, {1.0, 1.5,0.0,1.0} };
		double[] v1 = Rn.subtract(null, arc[1], arc[0]);
		double[] v2 = Rn.subtract(null, arc[2], arc[3]);
		FactoredMatrix M = new FactoredMatrix(P3.makeDirectIsometryFromFrames(null, arc[0], v1, new double[] {0,0,1,0}, arc[3], v2, new double[] {0,0,1,0}, Pn.EUCLIDEAN));
		double[] origin = Pn.dehomogenize(null, P2.pointFromLines(null, 
				P2.lineFromPoints(null, P2.projectP3ToP2(null, arc[0]), P2.projectP3ToP2(null, arc[1])),
				P2.lineFromPoints(null, P2.projectP3ToP2(null, arc[3]), P2.projectP3ToP2(null, arc[2]))
				));
		double[] t = M.getTranslation();
		double period = M.getRotationAngle();
		int n = 36;
		SceneGraphComponent e1 = new SceneGraphComponent();
		e1.setGeometry(IndexedFaceSetUtility.constructPolygon(conicalArc(origin, period, n, arc[0], arc[1])));
		SceneGraphComponent e2 = new SceneGraphComponent();
		e2.setGeometry(IndexedFaceSetUtility.constructPolygon(conicalArc(origin, period, n, arc[1], arc[2])));
		SceneGraphComponent e3 = new SceneGraphComponent();
		e3.setGeometry(IndexedFaceSetUtility.constructPolygon(conicalArc(origin, period, n, arc[2], arc[3])));
		SceneGraphComponent e4 = new SceneGraphComponent();
		e4.setGeometry(IndexedFaceSetUtility.constructPolygon(conicalArc(origin, period, n, arc[3], arc[0])));
		
		root.addChild(e1);
		root.addChild(e2);
		root.addChild(e3);
		root.addChild(e4);
		
		JRViewer.display(root);
	}
	
}
