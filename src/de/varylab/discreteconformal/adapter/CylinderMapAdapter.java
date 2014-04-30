package de.varylab.discreteconformal.adapter;

import static de.jreality.shader.CommonAttributes.EDGE_DRAW;
import static de.jreality.shader.CommonAttributes.VERTEX_DRAW;

import java.util.Comparator;
import java.util.Iterator;
import java.util.TreeSet;

import de.jreality.geometry.IndexedFaceSetUtility;
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
import de.varylab.discreteconformal.heds.CoVertex;

public class CylinderMapAdapter extends AbstractAdapter<double[]> {

	private Vertex<?,?,?> root = null;
	private double[] rootT = null;
	private int dir = 0;
	private double period = 1.0;
	
	public CylinderMapAdapter(Vertex<?, ?, ?> v, int direction, double p) {
		super(double[].class, true, false);
		root = v;
		this.dir = direction;
		period = p;
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
		return cylinder(period/(2*Math.PI), dir, vT);
	}

	private static double[] cylinder(double radius, int dir, double[] vT) {
		double 
			h = vT[(dir+1)%2], 
			phi = vT[dir];
		return new double[]{radius*Math.cos(phi/radius), radius*Math.sin(phi/radius), h};
	}

	public void initialize(CoVertex v1, int direction, double dist) {
		root = v1;
		dir = direction;
		period = dist;
		rootT = null;
	}
	
	public double[][] getArc(Edge<?,?,?> e, AdapterSet a, int n) {
		
		if(rootT == null) {
			rootT = a.getD(TexturePosition2d.class, root);
		}
		
		double[] start = Rn.subtract(null, a.getD(TexturePosition2d.class, e.getStartVertex()), rootT);
		double[] target = Rn.subtract(null, a.getD(TexturePosition2d.class, e.getTargetVertex()), rootT);
		
		double[][] curve = cylindricalArc(period/(2*Math.PI), dir, n, start, target);
		return curve;
	}

	private static class DirectionComparator implements Comparator<double[]> {

		private int dir = 0;
		
		public DirectionComparator(int d) {
			dir = d;
		}
		
		@Override
		public int compare(double[] o1, double[] o2) {
			return (int)Math.signum(o1[dir]-o2[dir]);
		}
		
	}
	
	private static double[][] cylindricalArc(double radius, int dir, int n, double[] start, double[] target) {
		if(start[dir] > target[dir]) {
			double[] tmp = start;
			start = target;
			target = tmp;
		}
		double[] dif = Rn.subtract(null, target, start);
		double d = 2*Math.PI*radius/n;
		double[] step = Rn.times(null, d/dif[dir], dif);
		
		
		
		int kstart = (int)(Math.ceil(start[dir]/d));
		int ktarget = (int)(Math.floor(target[dir]/d));
		TreeSet<double[]> domainVerts = new TreeSet<double[]>(new DirectionComparator(dir));
		domainVerts.add(start);
		domainVerts.add(target);
		
		for(int i = kstart; i <= ktarget; ++i) {
			double[] pt = Rn.add(null, start, Rn.times(null, i - start[dir]/d, step));
			domainVerts.add(pt);
		}
		double[][] curve = new double[domainVerts.size()][3];
		int i = 0;
		for (Iterator<double[]> iterator = domainVerts.iterator(); iterator.hasNext();) {
			curve[i++] = cylinder(radius,dir,iterator.next());
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
		double[][] arc = new double[][]{ {0.0,0.0},{Math.PI, 0.0}, {.5*Math.PI,.5} };
		
		int n = 36;
		double r = 0.8;
		SceneGraphComponent e1 = new SceneGraphComponent();
		e1.setGeometry(IndexedFaceSetUtility.constructPolygon(cylindricalArc(r, 0, n, arc[0], arc[1])));
		SceneGraphComponent e2 = new SceneGraphComponent();
		e2.setGeometry(IndexedFaceSetUtility.constructPolygon(cylindricalArc(r, 0, n, arc[1], arc[2])));
		SceneGraphComponent e3 = new SceneGraphComponent();
		e3.setGeometry(IndexedFaceSetUtility.constructPolygon(cylindricalArc(r, 0, n, arc[2], arc[0])));
		
		root.addChild(e1);
		root.addChild(e2);
		root.addChild(e3);
		
		JRViewer.display(root);
	}
}
