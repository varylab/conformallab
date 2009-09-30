package de.varylab.discreteconformal.plugin;

import static javax.swing.JOptionPane.WARNING_MESSAGE;
import geom3d.Point;
import geom3d.Triangle;

import java.awt.GridBagConstraints;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import de.jreality.plugin.basic.View;
import de.jreality.plugin.experimental.ManagedContent;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.jreality.adapter.Adapter;
import de.jtem.halfedgetools.plugin.HalfedgeConnectorPlugin;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.PluginInfo;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.varylab.discreteconformal.adapter.PositionAdapter;
import de.varylab.discreteconformal.adapter.TexCoordAdapter;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.bsp.HasPosition;
import de.varylab.discreteconformal.heds.bsp.KdTree;

public class SurfaceRemeshingPlugin extends ShrinkPanelPlugin implements ActionListener {

	public static abstract class MeshPattern {
		  

		
	}
	
	// plug-in connection
	private ManagedContent
		managedContent = null;
	private HalfedgeConnectorPlugin<CoVertex, CoEdge, CoFace, CoHDS>
		hcp = null;
	
	// ui components
	private GridBagConstraints
		gbc1 = new GridBagConstraints(),
		gbc2 = new GridBagConstraints();
	private JButton
		meshingButton = new JButton("Remesh Surface");
	
	public SurfaceRemeshingPlugin() {
		gbc1.insets = new Insets(2, 2, 2, 2);
		gbc1.gridwidth = GridBagConstraints.RELATIVE;
		gbc1.weightx = 0.0;
		gbc2.insets = new Insets(2, 2, 2, 2);
		gbc2.gridwidth = GridBagConstraints.REMAINDER;
		gbc2.weightx = 1.0;
		shrinkPanel.add(meshingButton);

		meshingButton.addActionListener(this);
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == meshingButton) {
			remeshSurface();
		}
	}
	
	
	private void remeshSurface() {
		Adapter a1 = new PositionAdapter();
		Adapter a2 = new TexCoordAdapter(false);
		CoHDS hds = new CoHDS();
		hcp.getHalfedgeContent(hds, a1, a2);
		if (hds == null) {
			return;
		}
		Window w = SwingUtilities.getWindowAncestor(shrinkPanel);
		if (hds.getVertex(0).getTextureCoord() == null) {
			JOptionPane.showMessageDialog(w, "Surface has no texture coordinates.", "Error", WARNING_MESSAGE);
			return;
		}
		
		double minX = hds.getVertex(0).getTextureCoord().x();
		double maxX = hds.getVertex(0).getTextureCoord().x();
		double minY = hds.getVertex(0).getTextureCoord().y();
		double maxY = hds.getVertex(0).getTextureCoord().y();
		for (CoVertex v : hds.getVertices()) {
			if (v.getTextureCoord().x() < minX) {
				minX = v.getTextureCoord().x(); 
			}
			if (v.getTextureCoord().x() > maxX) {
				maxX = v.getTextureCoord().x();
			}
			if (v.getTextureCoord().y() < minY) {
				minY = v.getTextureCoord().y(); 
			}
			if (v.getTextureCoord().y() > maxY) {
				maxY = v.getTextureCoord().y();
			}			
		}
		System.out.println("x[" + minX + ":" + maxX + "], y[" + minY + ":" + maxY + "]");
		double xSpan = maxX - minX;
		double ySpan = maxY - minY;
		int xRes = 100;
		int yRes = 100;
		double xStep = xSpan / xRes;
		double yStep = ySpan / yRes;
		CoHDS r = new CoHDS();
		for (int i = 0; i < xRes; i++) {
			for (int j = 0; j < yRes; j++) {
				CoVertex v = r.addNewVertex();
				double xPos = minX + i * xStep;
				double yPos = minY + j * yStep;
				Point pos = new Point(xPos,yPos,0);
				v.setPosition(pos);
			}
		}
		for (int i = 0; i < xRes - 1; i++) {
			for (int j = 0; j < yRes - 1; j++) {
				CoVertex v1 = r.getVertex(i*yRes + j); 
				CoVertex v2 = r.getVertex((i+1)*yRes + j); 
				CoVertex v3 = r.getVertex((i+1)*yRes + j + 1); 
				CoVertex v4 = r.getVertex(i*yRes + j + 1); 
				HalfEdgeUtils.constructFaceByVertices(r, v1, v2, v3, v4);
			}
		}
		PositionTexcoord[] vArr = new PositionTexcoord[hds.numVertices()];
		int index = 0; 
		for (CoVertex v : hds.getVertices()) {
			PositionTexcoord pt = new PositionTexcoord(v.getTextureCoord(), v.getPosition(), v);
			vArr[index++] = pt;
		}
		
		KdTree<PositionTexcoord> kdTree = new KdTree<PositionTexcoord>(vArr, 20, false);
		Set<CoVertex> cutted = new HashSet<CoVertex>(r.getVertices());
		for (CoVertex v : r.getVertices()) {
			for (PositionTexcoord ref : kdTree.collectKNearest(v, 15)) {
				Point patternPoint = v.getPosition();
				CoVertex refV = ref.getVertex();
				v.setTextureCoord(refV.getTextureCoord());
				for (CoFace f : HalfEdgeUtils.facesIncidentWithVertex(refV)) {
					List<CoVertex> b = HalfEdgeUtils.boundaryVertices(f);
					Triangle tTex = new Triangle(b.get(0).getTextureCoord(), b.get(1).getTextureCoord(), b.get(2).getTextureCoord());
					Triangle tPos = new Triangle(b.get(0).getPosition(), b.get(1).getPosition(), b.get(2).getPosition());
					if (isInTriangle(patternPoint, tTex, 0)) {
						Point bary = getBarycentic(patternPoint, tTex);
						Point newPos = getCoordinate(bary, tPos);
						v.setPosition(newPos);
						cutted.remove(v);
						break;
					}
				}
			}
		}
		
		
		hcp.updateHalfedgeContent(r, true, new PositionAdapter());
	}
	
	
	private boolean isInTriangle(Point p , Triangle t, double tol) {
		Point bary = getBarycentic(p, t);
		tol = -Math.abs(tol);
		return (bary.x() > tol && bary.y() > tol && bary.z() > tol);
	}
	
	
	/**
	 * Convert to barycentric
	 * @param p
	 * @param t
	 * @return
	 */
	private Point getBarycentic(Point p , Triangle t) {
		Point l = new Point();
		double x1 = t.getA().x();
		double y1 = t.getA().y();
		double x2 = t.getB().x();
		double y2 = t.getB().y();
		double x3 = t.getC().x();
		double y3 = t.getC().y();		
		double det = (x1 - x3)*(y2 - y3) - (y1 - y3)*(x2 - x3);
		l.setX(((y2 - y3)*(p.x() - x3) - (x2 - x3)*(p.y() - y3)) / det);
		l.setY(((x1 - x3)*(p.y() - y3) - (y1 - y3)*(p.x() - x3)) / det);
		l.setZ(1 - l.x() - l.y());
		return l;
	}
	
	private Point getCoordinate(Point b, Triangle t) {
		Point r = new Point();
		r.setX(b.x()*t.getA().x() + b.y()*t.getB().x() + b.z()*t.getC().x());
		r.setY(b.x()*t.getA().y() + b.y()*t.getB().y() + b.z()*t.getC().y());
		r.setZ(b.x()*t.getA().z() + b.y()*t.getB().z() + b.z()*t.getC().z());
		return r;
	}
	
	
	
	
	public class PositionTexcoord implements HasPosition {

		private Point
			t = null,
			p = null;
		private CoVertex	
			v = null;
		
		public PositionTexcoord(Point t, Point p, CoVertex v) {
			this.t = t;
			this.p = p;
			this.v = v;
		}
		
		@Override
		public Point getPosition() {
			return t;
		}

		@Override
		public void setPosition(Point p) {
			this.t.set(t.get());
		}
		
		public Point getPos() {
			return p;
		}
		
		public void setPos(Point t) {
			this.p.set(p.get());
		}
		
		public CoVertex getVertex() {
			return v;
		}
		
		public void setVertex(CoVertex vertex) {
			this.v = vertex;
		}
		
	}
	
	
	
	
	
	@SuppressWarnings("unchecked")
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hcp = c.getPlugin(HalfedgeConnectorPlugin.class);
		managedContent = c.getPlugin(ManagedContent.class);
	}
	
	
	@Override
	public void uninstall(Controller c) throws Exception {
		super.uninstall(c);
		managedContent.removeAll(getClass());
	}
	
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}

	@Override
	public PluginInfo getPluginInfo() {
		PluginInfo info = new PluginInfo("Surface Remeshing", "Stefan Sechelmann");
		return info;
	}

}
