package de.varylab.discreteconformal.plugin.schottky;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static javax.swing.JOptionPane.WARNING_MESSAGE;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

import no.uib.cipr.matrix.Vector;
import de.jreality.math.Rn;
import de.jreality.plugin.JRViewer;
import de.jreality.plugin.basic.View;
import de.jreality.ui.LayoutFactory;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.generic.Position3d;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeLayer;
import de.jtem.java2d.Viewer2D;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.mfc.field.Complex;
import de.jtem.mfc.group.Moebius;
import de.varylab.discreteconformal.functional.ConformalFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.StereographicProjectionAdapter;
import de.varylab.discreteconformal.math.ComplexUtility;
import de.varylab.discreteconformal.plugin.DiscreteConformalPlugin;
import de.varylab.discreteconformal.unwrapper.HyperbolicLayout;
import de.varylab.discreteconformal.unwrapper.HyperbolicUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.SphereUtility;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.NodeIndexComparator;
import de.varylab.discreteconformal.util.PathUtility;
import de.varylab.discreteconformal.util.Search;
import de.varylab.discreteconformal.util.SurgeryUtility;

public class SchottkyPlugin extends ShrinkPanelPlugin implements ActionListener {

	private HalfedgeInterface
		hif = null;
	private DiscreteConformalPlugin
		dcp = null;
	private SchottkyModeller
		schottkyModeller = new SchottkyModeller();
	private Viewer2D
		viewer = schottkyModeller.getViewer();
	private JButton
		generateButton = new JButton("Generate Surface"),
		toFuchsianButton = new JButton("Uniformize");
	private SpinnerNumberModel
		equalizationIterationsModel = new SpinnerNumberModel(10, 0, 100, 1),
		cirleResModel = new SpinnerNumberModel(20, 4, 1000, 1),
		extraPointsModel = new SpinnerNumberModel(400, 0, 10000, 1);
	private JSpinner
		equalizationIterationsSpinner = new JSpinner(equalizationIterationsModel),
		extraPointsSpinner = new JSpinner(extraPointsModel),
		circleResSpinner = new JSpinner(cirleResModel);
	private JCheckBox
		cutFromRootChecker = new JCheckBox("Use Cut Root");
	private Random
		rnd = new Random();
	private HalfedgeLayer
		schottkyLayer = null;
	
	
	public SchottkyPlugin() {
		GridBagConstraints c1 = LayoutFactory.createLeftConstraint();
		GridBagConstraints c2 = LayoutFactory.createRightConstraint();
		shrinkPanel.setTitle("Schottky Modeller");
		shrinkPanel.setFillSpace(true);
		viewer.setPreferredSize(new Dimension(240, 240));
		viewer.setMinimumSize(viewer.getPreferredSize());
		JScrollPane scrollProtector = new JScrollPane(viewer);
		scrollProtector.setPreferredSize(viewer.getPreferredSize());
		scrollProtector.setMinimumSize(scrollProtector.getPreferredSize());
		scrollProtector.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		scrollProtector.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_NEVER);
		
		shrinkPanel.setLayout(new GridBagLayout());
		c2.weighty = 1.0;
		shrinkPanel.add(scrollProtector, c2);
		c2.weighty = 0.0;
		shrinkPanel.add(new JLabel("Num Extra Points"), c1);
		shrinkPanel.add(extraPointsSpinner, c2);
		shrinkPanel.add(new JLabel("Point Equalizer Iterations"), c1);
		shrinkPanel.add(equalizationIterationsSpinner, c2);
		shrinkPanel.add(new JLabel("Circle Resolution"), c1);
		shrinkPanel.add(circleResSpinner, c2);
		shrinkPanel.add(cutFromRootChecker, c1);
		shrinkPanel.add(generateButton, c1);
		shrinkPanel.add(toFuchsianButton, c2);

		generateButton.addActionListener(this);
		toFuchsianButton.addActionListener(this);
		
		viewer.setScaleToolEnabled(true);
		viewer.setTranslateToolEnabled(true);
		viewer.setMenuToolEnabled(true);
	}
	
	
	@Override
	public void actionPerformed(ActionEvent event) {
		if (generateButton == event.getSource()) {
			List<SchottkyGenerator> pairs = schottkyModeller.getGenerators();
			Complex root = schottkyModeller.getBasePoint();
			Map<CoEdge, Double> lMap = new HashMap<CoEdge, Double>();
			Set<Set<CoEdge>> cylces = new HashSet<Set<CoEdge>>();
			CoHDS hds = new CoHDS();
			Map<CoVertex, Complex> mapCycleMap = new HashMap<CoVertex, Complex>();
			generate(hds, pairs, root, lMap, cylces, mapCycleMap);
			
			hif.set(hds);
			hif.addLayerAdapter(new SchottkyLengthAdapter(lMap), false);
		}
		if (toFuchsianButton == event.getSource()) {
			List<SchottkyGenerator> pairs = schottkyModeller.getGenerators();
			Complex root = schottkyModeller.getBasePoint();
			Map<CoEdge, Double> lMap = new HashMap<CoEdge, Double>();
			Set<Set<CoEdge>> cycles = new HashSet<Set<CoEdge>>();
			CoHDS hds = new CoHDS();
			Map<CoVertex, Complex> mapCycleMap = new HashMap<CoVertex, Complex>();
			CoVertex rootVertex = generate(hds, pairs, root, lMap, cycles, mapCycleMap);
			AdapterSet aSet = hif.getAdapters();
			aSet.add(new SchottkyLengthAdapter(lMap));
			
			int genus = HalfEdgeUtils.getGenus(hds);
			System.out.println("unwrapping surface of genus " + genus + "...");
			HyperbolicUnwrapperPETSc unwrapper = new HyperbolicUnwrapperPETSc();
			unwrapper.setCutAndLayout(false);
			unwrapper.setGradientTolerance(1E-8);
			unwrapper.setMaxIterations(2000);
			try {
				unwrapper.unwrap(hds, genus, aSet);
			} catch (Exception e1) {
				JOptionPane.showMessageDialog(shrinkPanel, e1.getMessage(), "Optimizer error", WARNING_MESSAGE);
				e1.printStackTrace();
				return;
			}

			CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
			for (Set<CoEdge> cycle : cycles) {
				Set<CoVertex> vCycle = PathUtility.getVerticesOnPath(cycle);
				CuttingUtility.cutAlongPath(cycle, cutInfo);
				for (CoVertex v : vCycle) {
					Complex z = mapCycleMap.get(v);
					CoVertex vv = cutInfo.vertexCopyMap.get(v);
					TopologyAlgorithms.removeVertex(vv);
//					aSet.set(Position.class, vv, new double[] {z.re, z.im});
				}
			}
			
			cutToSimplyConnected(hds, rootVertex, cutInfo);
			
			CoVertex layoutRoot = hds.getVertex(0);
			ConformalFunctional<CoVertex, CoEdge, CoFace> fun = unwrapper.getFunctional();
			Vector u = unwrapper.getUResult();
			layoutRoot = HyperbolicLayout.doLayout(hds, layoutRoot, fun,  u);
			
			dcp.createVisualization(hds, genus, cutInfo);
			dcp.updateSurface();
		}
	}
	
	
	private void cutToSimplyConnected(CoHDS hds, CoVertex cutRoot, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo) {
		List<List<CoEdge>> bc = HalfEdgeUtils.boundaryComponents(hds);
		if (cutRoot != null && cutFromRootChecker.isSelected()) {
			for (List<CoEdge> path : bc) {
				Set<CoVertex> cycle = PathUtility.getVerticesOnPath(path);
				List<CoEdge> cutPath = Search.getShortestPath(cutRoot, cycle, null);
				CuttingUtility.cutAlongPath(cutPath, cutInfo);
			}
		} else {
			while (bc.size() > 1) {
				Set<CoVertex> cycle1 = PathUtility.getVerticesOnPath(bc.get(0));
				List<CoEdge> path = null;
				Set<CoVertex> cycle2 = PathUtility.getVerticesOnPath(bc.get(1));
				
				Set<CoVertex> inter = new HashSet<CoVertex>(cycle1);
				inter.retainAll(cycle2);
				if (!inter.isEmpty()) {
					throw new RuntimeException("paths cannot intersect!");
				}
				for (CoVertex s : cycle1) {
					List<CoEdge> checkPath = Search.getShortestPath(s, cycle2, null);
					if (path == null) {
						path = checkPath;
					}
					if (checkPath.size() < path.size() && !checkPath.isEmpty()) {
						path = checkPath;
					}
				}
				if (path == null || path.isEmpty()) {
					throw new RuntimeException("no path between cycles found!");
				}
				CuttingUtility.cutAlongPath(path, cutInfo);
				bc = HalfEdgeUtils.boundaryComponents(hds);
			}
		}
	}
	
	
	
	private List<SchottkyCircle> getAllCircles(List<SchottkyGenerator> pairs) {
		List<SchottkyCircle> r = new LinkedList<SchottkyCircle>();
		for (SchottkyGenerator p : pairs) {
			SchottkyCircle sCircle = p.getCycle();
			SchottkyCircle tCircle = p.mapCircle(sCircle);
			r.add(sCircle);
			r.add(tCircle);
		}
		return r;
	}
	
	
	private void cutIdentificationHoles(CoHDS hds, List<ArrayList<CoVertex>> circles) {
		for (ArrayList<CoVertex> circle : circles) {
			int len = circle.size();
			for (int i = 0; i < len; i++) {
				CoVertex pre = circle.get((i + len - 1) % len);
				CoVertex act = circle.get((i + len + 0) % len);
				CoVertex pos = circle.get((i + len + 1) % len);
				List<CoEdge> star = incomingEdges(act);
				for (CoEdge e : star) {
					if (!e.isValid()) continue;
					if (e.getStartVertex() == pre) continue;
					if (e.getStartVertex() == pos) continue;
					if (circle.contains(e.getStartVertex())) {
						TopologyAlgorithms.removeEdge(e);
					}
				}
			}
		}
	}

	
	private CoVertex generate(CoHDS hds, List<SchottkyGenerator> pairs, Complex rootPos, Map<CoEdge, Double> lMap, Set<Set<CoEdge>> cyclesReturn, Map<CoVertex, Complex> mapCycleMap) {
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		int circleRes = cirleResModel.getNumber().intValue();
		int numExtraPoints = extraPointsModel.getNumber().intValue();
		
		Map<CoVertex, CoVertex> sMap = new HashMap<CoVertex, CoVertex>();
		Map<CoVertex, CoVertex> sInvMap = new HashMap<CoVertex, CoVertex>();
		Map<CoVertex, SchottkyGenerator> vertexPairMap = new HashMap<CoVertex, SchottkyGenerator>(); 
		Map<CoVertex, SchottkyGenerator> vertexPairInvMap = new HashMap<CoVertex, SchottkyGenerator>(); 
		Map<CoEdge, CoEdge> edgeMap = new HashMap<CoEdge, CoEdge>();
		Map<CoEdge, SchottkyGenerator> edgePairMap = new HashMap<CoEdge, SchottkyGenerator>();
		Map<CoEdge, SchottkyGenerator> edgePairInvMap = new HashMap<CoEdge, SchottkyGenerator>();
		
		
		// equalized extra vertices
		for (int i = 0; i < numExtraPoints; i++) {
			CoVertex v = hds.addNewVertex();
			double[] vPos = new double[] {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
			Rn.normalize(vPos, vPos);
			a.set(Position.class, v, vPos);
		}
		int numIterations = equalizationIterationsModel.getNumber().intValue();
		if (numIterations > 0) {
			Set<CoVertex> emptySet = Collections.emptySet();
			SphereUtility.equalizeSphereVertices(hds, emptySet, numIterations, 1E-6);
		}
		for (CoVertex v : hds.getVertices()) {
			double[] vPos = a.getD(Position.class, v);
			Complex z = ComplexUtility.stereographic(vPos);
			a.set(Position.class, v, new double[] {z.re, z.im, 0});
		}
		
//		// the root vertex 
		CoVertex root = hds.addNewVertex();
		a.set(Position.class, root, new double[]{rootPos.re, rootPos.im, 0});
		
		// add the vertices on the source and target circles
		Set<CoVertex> fixedVertices = new HashSet<CoVertex>();
		List<ArrayList<CoVertex>> sourceVertexCircles = new ArrayList<ArrayList<CoVertex>>();
		List<ArrayList<CoVertex>> targetVertexCircles = new ArrayList<ArrayList<CoVertex>>();
		for (SchottkyGenerator p : pairs) {
			ArrayList<CoVertex> sourceCircle = new ArrayList<CoVertex>();
			ArrayList<CoVertex> targetCircle = new ArrayList<CoVertex>();
			sourceVertexCircles.add(sourceCircle);
			targetVertexCircles.add(targetCircle);
			for (int i = 0; i < circleRes; i++) {
				SchottkyCircle c = p.getCycle();
				double phi = 2*i*PI / circleRes;
				double x = c.getRadius() * cos(phi) + c.getCenter().re;
				double y = c.getRadius() * sin(phi) + c.getCenter().im;
				Complex z = new Complex(x, y);
				Complex sz = p.getMoebius().applyTo(z);
				double[] zPos = new double[] {z.re, z.im, 0};
				double[] szPos = new double[] {sz.re, sz.im, 0};
				CoVertex v = hds.addNewVertex();
				CoVertex sv = hds.addNewVertex();
				sourceCircle.add(v);
				targetCircle.add(sv);
				a.set(Position.class, v, zPos);
				a.set(Position.class, sv, szPos);
				sMap.put(v, sv);
				sInvMap.put(sv, v);
				vertexPairMap.put(v, p);
				vertexPairInvMap.put(sv, p);
				fixedVertices.add(v);
				fixedVertices.add(sv);
			}
		}

		// exclude extra vertices from the circles
		for (CoVertex v : new HashSet<CoVertex>(hds.getVertices())) {
			if (sMap.containsKey(v) || sMap.containsValue(v)) {
				continue; // we are a circle vertex
			}
			for (SchottkyCircle c : getAllCircles(pairs)) {
				double[] vp = a.getD(Position3d.class, v); 
				Complex vz = new Complex(vp[0], vp[1]);
				if (c.isInside(vz, c.getRadius() * 1E-4)) {
					if (!v.isValid()) {
						continue;
					}
					hds.removeVertex(v);
				}
			}
		}
		
		// scale and project 
		Map<CoVertex, Complex> zMap = new HashMap<CoVertex, Complex>();
		for (CoVertex v : hds.getVertices()) {
			double[] p = a.getD(Position3d.class, v);
			Complex z = new Complex(p[0], p[1]);
			zMap.put(v, z); // store original positions
			double[] pProjected = ComplexUtility.inverseStereographic(z);
			a.set(Position.class, v, pProjected);
		}
		
		// create triangulation
		ConvexHull.convexHull(hds, a, 1E-8);
		List<ArrayList<CoVertex>> vertexCircles = new LinkedList<ArrayList<CoVertex>>();
		vertexCircles.addAll(sourceVertexCircles);
		vertexCircles.addAll(targetVertexCircles);
		cutIdentificationHoles(hds, vertexCircles);
		
		// back to C
		for (CoVertex v : hds.getVertices()) {
			Complex z = zMap.get(v);
			a.set(Position.class, v, new double[] {z.re, z.im});
		}
				
		schottkyLayer.addAdapter(new StereographicProjectionAdapter(), false);
		schottkyLayer.set(hds);
		
		// build edge identification opposites maps
		for (List<CoVertex> circle : sourceVertexCircles) {
			Set<CoVertex> vSet = new HashSet<CoVertex>(circle);
			for (CoVertex vt : circle) {
				CoVertex svt = sMap.get(vt);
				for (CoEdge e : incomingEdges(vt)) {
					CoVertex vs = e.getStartVertex();
					if (vSet.contains(vs)) {
						CoVertex svs = sMap.get(vs);
						for (CoEdge se : incomingEdges(svs)) {
							if (se.getStartVertex() == svt) {
								SchottkyGenerator p = vertexPairMap.get(vt);
								edgeMap.put(e, se);
								edgeMap.put(e.getOppositeEdge(), se.getOppositeEdge());
								edgeMap.put(se, e);
								edgeMap.put(se.getOppositeEdge(), e.getOppositeEdge());
								edgePairMap.put(e, p);
								edgePairMap.put(e.getOppositeEdge(), p);
								edgePairInvMap.put(se, p);
								edgePairInvMap.put(se.getOppositeEdge(), p);
							}
						}
					}
				}
			}
		}
		
		// we know how many edges there are on the circles
		assert edgeMap.size() == 2 * pairs.size() * circleRes * 2 : "lengths of cycles";
		
		// calculate length cross ratios
		Map<CoEdge, Double> crMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : hds.getEdges()) {
			CoVertex vi = e.getStartVertex();
			CoVertex vj = e.getTargetVertex();
			CoVertex vl = e.getNextEdge().getTargetVertex();
			CoVertex vk = e.getOppositeEdge().getNextEdge().getTargetVertex();
			Complex zi = zMap.get(vi);
			Complex zj = zMap.get(vj);
			Complex zl = zMap.get(vl);
			Complex zk = zMap.get(vk);
			
			Moebius pullTransform = null;
			// check if we have to pull a vertex position
			if (edgePairMap.containsKey(e)) {
				assert e.getLeftFace() == null || e.getRightFace() == null;
				// we are at a source circle edge
				SchottkyGenerator pair = edgePairMap.get(e);
				pullTransform = pair.getMoebius().invert();
			} 
			if (edgePairInvMap.containsKey(e)) {
				assert e.getLeftFace() == null || e.getRightFace() == null;
				// we are at a target circle edge
				SchottkyGenerator pair = edgePairInvMap.get(e);
				pullTransform = pair.getMoebius();
			}
			// pull one vertex from a different domain
			if (pullTransform != null) {
				CoEdge se = edgeMap.get(e);
				if (se.getLeftFace() == null) {
					CoVertex mappedL = se.getOppositeEdge().getNextEdge().getTargetVertex();
					zl = ComplexUtility.stereographic(a.getD(Position3d.class, mappedL));
					zl = pullTransform.applyTo(zl);
				} else {
					CoVertex mappedL = se.getNextEdge().getTargetVertex();
					zk = ComplexUtility.stereographic(a.getD(Position3d.class, mappedL));
					zk = pullTransform.applyTo(zk);
				}
			}
			double l_ik = zi.dist(zk); 
			double l_kj = zk.dist(zj); 
			double l_jl = zj.dist(zl); 
			double l_li = zl.dist(zi); 
			double cr = (l_kj * l_li) / (l_ik * l_jl);
			// store the cross ratio
			crMap.put(e, cr);
		}
		
		// store to be lost coordinates
		for (CoEdge e : edgeMap.keySet()) {
			CoEdge se = edgeMap.get(e);
//			Complex mape = zMap.get(e.getStartVertex());
			Complex mapse = zMap.get(se.getStartVertex());
			mapCycleMap.put(e.getStartVertex(), mapse);
//			mapCycleMap.put(se.getStartVertex(), mapse);
		}

		// identify circles
		try {
			Set<CoEdge> mappedEdgesSet = new TreeSet<CoEdge>(new NodeIndexComparator<CoEdge>());
			mappedEdgesSet.addAll(edgeMap.keySet());
			for (CoEdge e : mappedEdgesSet) {
				if (!e.isValid() || e.getLeftFace() != null) {
					continue;
				}
				CoEdge se = edgeMap.get(e);
				try {
					System.out.println("Identify " + e + " " + se);
					Set<CoEdge> cycle = SurgeryUtility.glueAlongBoundaries(e, se);
					cyclesReturn.add(cycle);
				} catch (RuntimeException re) {
					re.printStackTrace();
				}
			}
		} catch (Exception e) {
			System.err.println("Invalid Schottky configuration.");
		}
		
		// check if identification was correct
		for (CoEdge e : edgeMap.keySet()) {
			if (!e.isValid()) continue;
			CoEdge se = edgeMap.get(e);
			assert e == se.getOppositeEdge() : "opposite";
		}
		// assert mapped and opposite cross ratios are the same
		for (CoEdge e : edgeMap.keySet()) {
			CoEdge se = edgeMap.get(e);
			double cr = crMap.get(e);
			double scr = crMap.get(se);
			assert abs(cr - scr) < 1E-8 : "Mapped cross ratio: " + cr + " - " + scr;
		}
		for (CoEdge e : hds.getEdges()) {
			CoEdge opp = e.getOppositeEdge();
			double cr = crMap.get(e);
			double scr = crMap.get(opp);
			assert abs(cr - scr) < 1E-8 : "Opposite cross ratios: " + cr + " - " + scr;
		}
		// check cross-ratio product around each vertex
		for (CoVertex v : hds.getVertices()) {
			double pr = 1.0;
			for (CoEdge e : incomingEdges(v)) {
				pr *= crMap.get(e);
			}
			assert abs(pr - 1.0) < 1E-8 : "Cross-ratio product around vertex not 1.0";
		}
		
		// define lengths from cross-ratios
		Map<CoEdge, Double> aMap = new HashMap<CoEdge, Double>();
		for (CoVertex v : hds.getVertices()) {
			double ai = 100.0;
			for (CoEdge e : incomingEdges(v)) {
				double q = crMap.get(e);
				ai *= q;
				aMap.put(e, ai);
			}
		}
		for (CoEdge e : hds.getEdges()) {
			double ai = aMap.get(e);
			double aj = aMap.get(e.getPreviousEdge());
			double l = sqrt(1 / (ai * aj));
			lMap.put(e, l);
		}
		
		// check opposite lengths
		for (CoEdge e : hds.getEdges()) {
			double l = lMap.get(e);
			double lo = lMap.get(e.getOppositeEdge());
			assert abs(l - lo) < 1E-8 : "lengths";
		}

		System.out.println("Generated surface of genus " + HalfEdgeUtils.getGenus(hds));
		return root.isValid() ? root : null;
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
		dcp = c.getPlugin(DiscreteConformalPlugin.class);
		schottkyLayer = new HalfedgeLayer(hif);
		schottkyLayer.setName("Schottky Data");
		hif.addLayer(schottkyLayer);
	}
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}
	
	
	public static void main(String[] args) {
		NativePathUtility.set("native");
		JRViewer v = new JRViewer();
		v.setPropertiesFile("Schottky.jrw");
		v.setPropertiesResource(SchottkyPlugin.class, "Schottky.jrw");
		v.addContentUI();
		v.addBasicUI();
		v.registerPlugin(SchottkyPlugin.class);
		v.startup();
	}
	
}
