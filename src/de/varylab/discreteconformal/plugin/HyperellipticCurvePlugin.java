package de.varylab.discreteconformal.plugin;

import static de.jreality.ui.LayoutFactory.createLeftConstraint;
import static de.jreality.ui.LayoutFactory.createRightConstraint;
import static de.varylab.discreteconformal.math.ComplexUtility.inverseStereographic;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

import de.jreality.math.Pn;
import de.jreality.plugin.basic.View;
import de.jtem.blas.ComplexMatrix;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.HalfedgeSelection;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.PluginInfo;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.surface.BranchPoint;
import de.jtem.riemann.theta.SiegelReduction;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.EllipticImageGenerator.PathVisualizer;
import de.varylab.discreteconformal.plugin.hyperelliptic.Curve;
import de.varylab.discreteconformal.plugin.hyperelliptic.CurveChangeEvent;
import de.varylab.discreteconformal.plugin.hyperelliptic.CurveChangeEvent.EventType;
import de.varylab.discreteconformal.plugin.hyperelliptic.CurveChangeListener;
import de.varylab.discreteconformal.plugin.hyperelliptic.CurveEditor;
import de.varylab.discreteconformal.plugin.image.ImageHook;
import de.varylab.discreteconformal.unwrapper.SphereUtility;
import de.varylab.discreteconformal.unwrapper.SphericalNormalizerPETc;
import de.varylab.discreteconformal.util.HyperellipticUtility;
import de.varylab.discreteconformal.util.SimpleMatrixPrintUtility;
import de.varylab.mtjoptimization.NotConvergentException;

public class HyperellipticCurvePlugin extends ShrinkPanelPlugin implements
		CurveChangeListener, ActionListener {

	private HalfedgeInterface
		hif = null;
	private Random
		rnd = new Random();
	
	private CurveEditor editor;
	
	private JPanel
		geometryPanel = new JPanel();
	private SpinnerNumberModel
		extraPointsModel = new SpinnerNumberModel(100, 0, 10000, 1),
		equalizationIterationsModel = new SpinnerNumberModel(10, 0, 100, 1);
	private JSpinner
		equalizationIterationsSpinner = new JSpinner(equalizationIterationsModel),
		extraPointsSpinner = new JSpinner(extraPointsModel);
	private JCheckBox
		showcutPathsChecker = new JCheckBox("Show Cut Paths"),
		normalizerBranchPointPositionsChecker = new JCheckBox("Normalize Branch Points", true);
	private JButton
		createButton = new JButton("Create Triangulated Surface");
	
	public CurveEditor getEditor() {
		return editor;
	}

	public HyperellipticCurvePlugin() {
		Curve c = new Curve(1);
		c.update();
		setCurve(c);

		editor.setPreferredSize(new Dimension(240, 240));
		editor.setMinimumSize(editor.getPreferredSize());
		JScrollPane scrollProtector = new JScrollPane(editor);
		scrollProtector.setPreferredSize(editor.getPreferredSize());
		scrollProtector.setMinimumSize(scrollProtector.getPreferredSize());
		scrollProtector.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		scrollProtector.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_NEVER);
		setInitialPosition(SHRINKER_RIGHT);
		GridBagConstraints c1 = createLeftConstraint();
		GridBagConstraints c2 = createRightConstraint();
		shrinkPanel.setTitle("Hyperelliptic Curve Plugin");
		c2.weighty = 1.0;
		shrinkPanel.add(scrollProtector, c2);
		
		geometryPanel.setBorder(BorderFactory.createTitledBorder("Triangulated Surface"));
		geometryPanel.setLayout(new GridBagLayout());
		geometryPanel.add(new JLabel("Extra Random Points"), c1);
		geometryPanel.add(extraPointsSpinner, c2);
		geometryPanel.add(new JLabel("Point Equalizer Iterations"), c1);
		geometryPanel.add(equalizationIterationsSpinner, c2);
		geometryPanel.add(showcutPathsChecker, c2);
		geometryPanel.add(normalizerBranchPointPositionsChecker, c2);
		geometryPanel.add(createButton, c2);
		shrinkPanel.add(geometryPanel, c2);
		
		createButton.addActionListener(this);
	}

	public Curve getCurve() {
		return editor.getCurve();
	}

	public void setCurve(Curve curve) {

		curve.setEps(1E-15);
		
		if (editor != null) {
			curve.setCurveChangeListeners(getCurve().getCurveChangeListeners());
			editor.setCurve(curve);
		} else {
			this.editor = new CurveEditor(curve);
			getCurve().addCurveChangeListener(this);
		}
		update();

		CurveChangeEvent e = new CurveChangeEvent(getCurve(), this,
				EventType.NEW_CURVE_SET);
		getCurve().fireCurveChangeEvent(e);

	}

	@Override
	public void actionPerformed(ActionEvent e) {
		AdapterSet a = hif.getAdapters();
		// create a triangulated surface for the active curve
		if (createButton == e.getSource()) {
			CoHDS hds = new CoHDS();
			Set<CoVertex> branchVertices = new HashSet<CoVertex>();
			for (Complex z : getBranchPoints()) {
				double[] p = inverseStereographic(z);
				CoVertex v = hds.addNewVertex();
				a.set(Position.class, v, p);
				branchVertices.add(v);
			}
			// add north pole if needed
			if (branchVertices.size() % 2 == 1) {
				CoVertex v = hds.addNewVertex();
				a.set(Position.class, v, new double[] {0, 0, 1});
				branchVertices.add(v);
			}
			try {
				if (normalizerBranchPointPositionsChecker.isSelected()) {
					for (CoVertex v : hds.getVertices()) v.T = v.P.clone();
					SphericalNormalizerPETc.normalize(hds);
					for (CoVertex v : hds.getVertices()) v.P = v.T.clone();
				}
			} catch (NotConvergentException e1) {
				System.err.println("could nor normalize branch points " + e1.getLocalizedMessage());
			}
			int numextra = extraPointsModel.getNumber().intValue();
			List<CoVertex> extraVertices = new LinkedList<CoVertex>();
			// additional points
			for (int j = 0; j < numextra; j++) {
				CoVertex v = hds.addNewVertex();
				double[] p = new double[] {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian(), 1.0};
				Pn.setToLength(p, p, 1, Pn.EUCLIDEAN);
				a.set(Position.class, v, p);
				extraVertices.add(v);
			}
			// equalization
			int numIterations = equalizationIterationsModel.getNumber().intValue();
			if (numIterations > 0) {
				SphereUtility.equalizeSphereVertices(hds, branchVertices, numIterations, 1E-6);
			}
			Set<CoEdge> glueSet = new HashSet<CoEdge>();
			int[] branchIndices = new int[branchVertices.size()];
			int i = 0;
			for (CoVertex bv : branchVertices) {
				branchIndices[i++] = bv.getIndex();
			}
			HyperellipticUtility.generateHyperellipticImage(hds, 0, glueSet, branchIndices);
			if (showcutPathsChecker.isSelected()) {
				PathVisualizer pathVisualizer = new PathVisualizer();
				for (CoEdge pe : glueSet) {
					pathVisualizer.add(pe);
					pathVisualizer.add(pe.getOppositeEdge());
				}
				// show the result
				hif.addLayerAdapter(pathVisualizer, true);
			}
			hif.set(hds);
			HalfedgeSelection branchSelection = new HalfedgeSelection(branchVertices);
			hif.setSelection(branchSelection);
		}
	}
	
	@Override
	public void curveChanged(CurveChangeEvent e) {
		if (e.type == EventType.CURVE_CHANGED)
			updatePeriodMatrix();
	}

	private ComplexMatrix getNormalizedPeriodMatrix() {
		final ComplexMatrix PeriodMatrix = editor.getCurve().getPeriodMatrix().copy();
		SiegelReduction siegel = new SiegelReduction(PeriodMatrix);
		return siegel.getReducedPeriodMatrix();
	}

	public void update() {
		editor.getCurve().setEps(1E-16);
		editor.update();
	}

	private void updatePeriodMatrix() {
		System.out.println(SimpleMatrixPrintUtility.toString(getNormalizedPeriodMatrix(), 20));
	}

	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
	}
	
	@Override
	public PluginInfo getPluginInfo() {
		PluginInfo info = super.getPluginInfo();
		info.name = "Curve";
		info.vendorName = "";
		info.isDynamic = true;
		info.icon = ImageHook.getIcon("wenteTorus16.png");
		return info;
	}

	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}
	
	public Complex[] getBranchPoints(){
		BranchPoint[] bp= editor.getCurve().getBranchPoints();
		Complex[] copy= new Complex[bp.length];
		for (int i = 0; i < copy.length; i++) {
			copy[i]= bp[i].getCoords();
		}
		return copy;
	}

}