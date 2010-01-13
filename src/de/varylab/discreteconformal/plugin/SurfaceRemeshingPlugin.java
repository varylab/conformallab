package de.varylab.discreteconformal.plugin;

import static javax.swing.JOptionPane.WARNING_MESSAGE;
import geom3d.Point;
import geom3d.Triangle;

import java.awt.GridBagConstraints;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;

import de.jreality.math.FactoredMatrix;
import de.jreality.math.Matrix;
import de.jreality.math.Pn;
import de.jreality.plugin.basic.View;
import de.jreality.plugin.content.ContentAppearance;
import de.jreality.plugin.experimental.ManagedContent;
import de.jreality.ui.AppearanceInspector;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.jreality.adapter.Adapter;
import de.jtem.halfedgetools.plugin.HalfedgeConnectorPlugin;
import de.jtem.halfedgetools.util.HalfEdgeUtilsExtra;
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

	private enum Pattern {
		Quads,
		Hexagons;
	}
	
	private ContentAppearance
		contentAppearance = null;
	private JComboBox
		patternCombo = new JComboBox(Pattern.values());
	private SpinnerNumberModel
		lookUpHeuristicModel = new SpinnerNumberModel(3, 1, 100, 1),
		resolutionModel = new SpinnerNumberModel(50, 1, 10000000, 1);
	private JSpinner
		lookUpHeuristicSpinner = new JSpinner(lookUpHeuristicModel),
		resolutionSpinner = new JSpinner(resolutionModel);
	
	// plug-in connection
	private ManagedContent
		managedContent = null;
	private HalfedgeConnectorPlugin
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
		gbc1.fill = GridBagConstraints.BOTH;
		gbc2.insets = new Insets(2, 2, 2, 2);
		gbc2.gridwidth = GridBagConstraints.REMAINDER;
		gbc2.weightx = 1.0;
		gbc2.fill = GridBagConstraints.BOTH;

		shrinkPanel.add(new JLabel("Pattern"), gbc1);
		shrinkPanel.add(patternCombo, gbc2);
		shrinkPanel.add(new JLabel("Lookup Heuristic"), gbc1);
		shrinkPanel.add(lookUpHeuristicSpinner, gbc2);
		shrinkPanel.add(new JLabel("Resolution"), gbc1);
		shrinkPanel.add(resolutionSpinner, gbc2);
		
		shrinkPanel.add(meshingButton, gbc2);

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
		
		AppearanceInspector ai = contentAppearance.getAppearanceInspector();
		Matrix texMatrix = ai.getTextureMatrix();
		Matrix texInvMatrix = texMatrix.getInverse();
		FactoredMatrix fm = new FactoredMatrix(texMatrix, Pn.EUCLIDEAN);
		double[] translation = fm.getTranslation();
		
		for (CoVertex v : hds.getVertices()) {
			v.setTextureCoord(transformTexCoord(v.getTextureCoord(), texMatrix));
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

		minX = Math.floor(minX) - 1;
		maxX = Math.ceil(maxX) + 1;
		minY = Math.floor(minY) - 1;
		maxY = Math.ceil(maxY) + 1;
		
		double xSpan = maxX - minX;
		double ySpan = maxY - minY;

		CoHDS r = new CoHDS();
		// create pattern
		switch (getPattern()) {
		case Quads:
			int xRes = (int)Math.ceil(xSpan * 2);
			double xStep = 0.5;
			int yRes = (int)Math.ceil(ySpan * 2); 
			double yStep = 0.5;
			double xOffset = translation[0] % 1.0;
			double yOffset = translation[1] % 1.0;
			for (int i = 0; i < xRes; i++) {
				for (int j = 0; j < yRes; j++) {
					CoVertex v = r.addNewVertex();
					double xPos = minX + i * xStep;
					double yPos = minY + j * yStep;
					Point pos = new Point(xPos, yPos, 1);
					v.setPosition(pos);
					Point tex = new Point(xPos, yPos, 1);
					v.setTextureCoord(transformTexCoord(tex, texInvMatrix));
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
			break;
		case Hexagons:
			xRes = (int)Math.ceil(xSpan * 6/3.0);
			xStep = 3/6.0;
			yRes = (int)Math.ceil(ySpan * 2);
			yStep = 0.5;
			double radius = 1.0/3.0;
			xOffset = translation[0] + 0.75*radius % 1.0;
			yOffset = translation[1] % 1.0;
			for (int i = 0; i < xRes; i++) {
				for (int j = 0; j < yRes; j++) {
					CoVertex v = r.addNewVertex();
					double xPos = minX + i * xStep + xOffset;
					double yPos = minY + j * yStep + yOffset;
					double move = (j % 2 == 0) ? 1 : -1;
					move *= (i % 2 == 0) ? 1 : -1;
					xPos += move*radius/4;
					Point pos = new Point(xPos, yPos, 1);
					v.setPosition(pos);
					Point tex = new Point(xPos, yPos, 1);
					v.setTextureCoord(transformTexCoord(tex, texInvMatrix));
				}
			}
			
			for (int i = 0; i < xRes - 1; i++) {
				for (int j = 0; j < yRes - 2; j+=2) {
					int colStep = i % 2 == 0 ? 0 : 1;
					if (j + 2 + colStep >= yRes) {
						continue;
					}
					CoVertex v1 = r.getVertex(i*yRes + j + colStep); 
					CoVertex v2 = r.getVertex((i+1)*yRes + j + colStep); 
					CoVertex v3 = r.getVertex((i+1)*yRes + j + 1 + colStep); 
					CoVertex v4 = r.getVertex((i+1)*yRes + j + 2 + colStep); 
					CoVertex v5 = r.getVertex(i*yRes + j + 2 + colStep); 
					CoVertex v6 = r.getVertex(i*yRes + j + 1 + colStep); 
					HalfEdgeUtils.constructFaceByVertices(r, v1, v2, v3, v4, v5, v6);
				}
			}
		}
		
		
		
		PositionTexcoord[] vArr = new PositionTexcoord[hds.numVertices()];
		int index = 0; 
		for (CoVertex v : hds.getVertices()) {
			PositionTexcoord pt = new PositionTexcoord(v.getTextureCoord(), v.getPosition(), v);
			vArr[index++] = pt;
		}
		
		int lookUp = lookUpHeuristicModel.getNumber().intValue();
		Map<CoVertex, CoFace> texFaceMap = new HashMap<CoVertex, CoFace>();
		KdTree<PositionTexcoord> kdTree = new KdTree<PositionTexcoord>(vArr, 20, false);
		Set<CoVertex> cutted = new HashSet<CoVertex>(r.getVertices());
		for (CoVertex v : r.getVertices()) {
			Point patternPoint = v.getPosition();
			Vector<PositionTexcoord> posTexVec = kdTree.collectKNearest(v, lookUp);
			Set<CoFace> checkFaces = new HashSet<CoFace>();
			for (PositionTexcoord posTex : posTexVec) {
				checkFaces.addAll(HalfEdgeUtils.facesIncidentWithVertex(posTex.getVertex()));
			}
			for (CoFace f : checkFaces) {
				List<CoVertex> b = HalfEdgeUtils.boundaryVertices(f);
				Triangle tTex = new Triangle(b.get(0).getTextureCoord(), b.get(1).getTextureCoord(), b.get(2).getTextureCoord());
				Triangle tPos = new Triangle(b.get(0).getPosition(), b.get(1).getPosition(), b.get(2).getPosition());
				if (isInTriangle(patternPoint, tTex, 0)) {
					Point bary = getBarycentic(patternPoint, tTex);
					Point newPos = getCoordinate(bary, tPos);
					v.setPosition(newPos);
					cutted.remove(v);
					texFaceMap.put(v, f);
					break;
				}
			}
		}
		
		
		Set<CoVertex> overlap = new HashSet<CoVertex>();
		for (CoVertex v : r.getVertices()) {
			if (cutted.contains(v)) {
				continue;
			}
			CoFace mapFace = texFaceMap.get(v);
			List<CoEdge> star = HalfEdgeUtilsExtra.getEdgeStar(v);
			for (CoEdge e : star) {
				CoVertex incident = e.getStartVertex();
				if (cutted.contains(incident)) {
					CoFace overlapFace = e.getLeftFace();
					if (overlapFace == null) {
						continue;
					}
					List<CoVertex> b = HalfEdgeUtils.boundaryVertices(overlapFace);
					for (CoVertex bv : b) {
						if (cutted.contains(bv)) {
							overlap.add(bv);
							texFaceMap.put(bv, mapFace);
						}
					}
				}
			}
		}
		cutted.removeAll(overlap);
		for (CoVertex v : overlap) {
			CoFace f = texFaceMap.get(v);
			List<CoVertex> b = HalfEdgeUtils.boundaryVertices(f);
			Triangle tTex = new Triangle(b.get(0).getTextureCoord(), b.get(1).getTextureCoord(), b.get(2).getTextureCoord());
			Triangle tPos = new Triangle(b.get(0).getPosition(), b.get(1).getPosition(), b.get(2).getPosition());
			Point bary = getBarycentic(v.getPosition(), tTex);
			Point newPos = getCoordinate(bary, tPos);
			v.setPosition(newPos);
		}
		
		
		Set<CoFace> deleteFaces = new HashSet<CoFace>();
		Set<CoEdge> deleteEdges = new HashSet<CoEdge>();
		for (CoVertex v : cutted) {
			for (CoEdge e : HalfEdgeUtils.incomingEdges(v)) {
				deleteEdges.add(e);
				deleteEdges.add(e.getOppositeEdge());
			}
			deleteFaces.addAll(HalfEdgeUtils.facesIncidentWithVertex(v));
		}
		for (CoVertex v : cutted) {
			r.removeVertex(v);
		}
		for (CoEdge e : deleteEdges) {
			r.removeEdge(e);
		}
		for (CoFace f : deleteFaces) {
			r.removeFace(f);
		}
		hcp.updateHalfedgeContent(r, true, new PositionAdapter(), new TexCoordAdapter(true));
	}
	
	
	
	private Point transformTexCoord(Point tex, Matrix M) {
		double[] tmp = {0,0,0,1};
		tmp[0] = tex.z() != 0 ? tex.x() / tex.z() : tex.x();
		tmp[1] = tex.z() != 0 ? tex.y() / tex.z() : tex.y();
		tmp[2] = 0;
		tmp[3] = tex.z() != 0 ? tex.z() : 1;
		double[] texCoord = M.multiplyVector(tmp);
		return new Point(texCoord[0] / texCoord[3], texCoord[1] / texCoord[3], 1.0);
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
	
	
	
	
	public class BaryFace implements HasPosition {
		
		private Point
			bary = null;
		private CoFace
			face = null;
		
		public BaryFace(CoFace f) {
			this.face = f;
			List<CoVertex> b = HalfEdgeUtils.boundaryVertices(f);
			Triangle t = new Triangle(b.get(0).getTextureCoord(), b.get(1).getTextureCoord(), b.get(2).getTextureCoord());
			bary = t.getBaryCenter();
		}
		
		
		@Override
		public Point getPosition() {
			return bary;
		}
		@Override
		public void setPosition(Point p) {
			bary.set(p.get());
		}

		public CoFace getFace() {
			return face;
		}
		
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
	
	
	
	private Pattern getPattern() {
		return (Pattern)patternCombo.getSelectedItem();
	}
	
	
	
	
	@Override 
	public void install(Controller c) throws Exception {
		super.install(c);
		hcp = c.getPlugin(HalfedgeConnectorPlugin.class);
		managedContent = c.getPlugin(ManagedContent.class);
		contentAppearance = c.getPlugin(ContentAppearance.class);
	}
	
	
	@Override
	public void uninstall(Controller c) throws Exception {
		super.uninstall(c);
		managedContent.removeAll(getClass());
	} 
	
	
	@Override
	public void storeStates(Controller c) throws Exception {
		super.storeStates(c);
		c.storeProperty(getClass(), "resolution", resolutionModel.getNumber().intValue());
		c.storeProperty(getClass(), "lookUpHeuristic", lookUpHeuristicModel.getNumber().intValue());
		c.storeProperty(getClass(), "pattern", patternCombo.getSelectedIndex());
	}
	
	@Override
	public void restoreStates(Controller c) throws Exception {
		super.restoreStates(c);
		resolutionModel.setValue(c.getProperty(getClass(), "resolution", resolutionModel.getNumber().intValue()));
		lookUpHeuristicModel.setValue(c.getProperty(getClass(), "lookUpHeuristic", lookUpHeuristicModel.getNumber().intValue()));
		patternCombo.setSelectedIndex(c.getProperty(getClass(), "pattern", 0));
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
