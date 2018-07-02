package de.varylab.discreteconformal.plugin;

import static de.jreality.ui.LayoutFactory.createLeftConstraint;
import static de.jreality.ui.LayoutFactory.createRightConstraint;
import static de.varylab.discreteconformal.math.ComplexUtility.inverseStereographic;
import static de.varylab.discreteconformal.math.ComplexUtility.stereographic;
import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

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
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

import de.jreality.math.Matrix;
import de.jreality.math.MatrixBuilder;
import de.jreality.math.Pn;
import de.jreality.plugin.basic.View;
import de.jtem.blas.ComplexMatrix;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.selection.Selection;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.PluginInfo;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.surface.BranchPoint;
import de.jtem.riemann.theta.SiegelReduction;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.HyperEllipticAlgebraicCurve;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.hyperelliptic.Curve;
import de.varylab.discreteconformal.plugin.hyperelliptic.CurveChangeEvent;
import de.varylab.discreteconformal.plugin.hyperelliptic.CurveChangeEvent.EventType;
import de.varylab.discreteconformal.plugin.hyperelliptic.CurveChangeListener;
import de.varylab.discreteconformal.plugin.hyperelliptic.CurveEditor;
import de.varylab.discreteconformal.plugin.image.ImageHook;
import de.varylab.discreteconformal.unwrapper.SphereUtility;
import de.varylab.discreteconformal.unwrapper.SphericalNormalizerPETSc;
import de.varylab.discreteconformal.util.HyperellipticUtility;
import de.varylab.discreteconformal.util.SimpleMatrixPrintUtility;

public class HyperellipticCurvePlugin extends ShrinkPanelPlugin implements CurveChangeListener, ActionListener {

	private static Logger
		log = Logger.getLogger(HyperellipticCurvePlugin.class.getName());
	private HalfedgeInterface
		hif = null;
	private ConformalDataPlugin
		conformalDataPlugin = null;
	private Random
		rnd = new Random();
	
	private CurveEditor editor;
	
	private JPanel
		geometryPanel = new JPanel();
	private SpinnerNumberModel
		randomSeedModel = new SpinnerNumberModel(0, 0, 10000000, 1),
		extraPointsModel = new SpinnerNumberModel(0, 0, 10000, 1),
		extraPointsAtBranchModel = new SpinnerNumberModel(100, 0, 10000, 1),
		equalizationIterationsModel = new SpinnerNumberModel(0, 0, 1000, 1);
	private JSpinner
		randomSeedSpinner = new JSpinner(randomSeedModel),
		equalizationIterationsSpinner = new JSpinner(equalizationIterationsModel),
		extraPointsSpinner = new JSpinner(extraPointsModel),
		extraPointsAtBranchSpinner = new JSpinner(extraPointsAtBranchModel);
	private JCheckBox
		fibonacciChecker = new JCheckBox("Use Fibonacci Points", true);
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
		
		c2.weighty = 0.0;
		geometryPanel.setBorder(BorderFactory.createTitledBorder("Triangulated Surface"));
		geometryPanel.setLayout(new GridBagLayout());
		geometryPanel.add(new JLabel("Random Seed"), c1);
		geometryPanel.add(randomSeedSpinner, c2);
		geometryPanel.add(new JLabel("Extra Random Points"), c1);
		geometryPanel.add(extraPointsSpinner, c2);
		geometryPanel.add(new JLabel("Extra Random Points At Branches"), c1);
		geometryPanel.add(extraPointsAtBranchSpinner, c2);
		geometryPanel.add(fibonacciChecker, c2);
		geometryPanel.add(new JLabel("Point Equalizer Iterations"), c1);
		geometryPanel.add(equalizationIterationsSpinner, c2);
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
			update();
		} else {
			this.editor = new CurveEditor(curve);
			update();
			getCurve().addCurveChangeListener(this);
		}

		CurveChangeEvent e = new CurveChangeEvent(getCurve(), this,
				EventType.NEW_CURVE_SET);
		getCurve().fireCurveChangeEvent(e);

	}

	@Override
	public void actionPerformed(ActionEvent e) {
		rnd.setSeed(randomSeedModel.getNumber().intValue());
		if (createButton == e.getSource()) {
			Set<CoVertex> branchVertices = new HashSet<>();
			CoHDS hds = HyperellipticCurvePlugin.generateCurve(
				getBranchPoints(),
				extraPointsModel.getNumber().intValue(),
				extraPointsAtBranchModel.getNumber().intValue(),
				fibonacciChecker.isSelected(),
				equalizationIterationsModel.getNumber().intValue(),
				rnd,
				hif.getAdapters(),
				branchVertices
			);
			
			int genus = HalfEdgeUtils.getGenus(hds);
			HyperEllipticAlgebraicCurve heac = DataUtility.toHyperEllipticAlgebraicCurve("Hyperelliptic Curve g" + genus, getCurve());
			conformalDataPlugin.addData(heac);
			hif.set(hds);
			Selection branchSelection = new Selection(branchVertices);
			hif.setSelection(branchSelection);
		}
	}
	
	public static CoHDS generateCurve(
		Complex[] branchPoints,
		int numextra,
		int numextrabranch,
		boolean fibonacci,
		int numEqualizerIterations,
		Random rnd,
		AdapterSet a,
		Set<CoVertex> branchVerticesOUT
	) {
		CoHDS hds = new CoHDS();
		Set<CoVertex> branchVertices = new HashSet<>();
		for (Complex z : branchPoints) {
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
			SphericalNormalizerPETSc.normalize(hds, a, Position4d.class, Position.class);
		} catch (Exception e1) {
			log.warning("could not normalize branch points " + e1.getLocalizedMessage());
		}
		
		List<CoVertex> extraVertices = new LinkedList<CoVertex>();
		// additional points
		for (int j = 0; j < numextra; j++) {
			CoVertex v = hds.addNewVertex();
			double[] p;
			if (fibonacci) {
				double phi = acos(1 - 2 * (j + 0.5) / numextra);
				double theta = PI * (1 + sqrt(5)) * j;
				p = new double[] {cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi), 1.0};
			} else {
				p = new double[] {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian(), 1.0};
				Pn.setToLength(p, p, 1, Pn.EUCLIDEAN);
			}
			a.set(Position.class, v, p);
			extraVertices.add(v);
		}
		
		for (CoVertex bv : branchVertices) {
			double[] branchPos = a.getD(Position4d.class, bv); 
			Matrix T = MatrixBuilder.euclidean().rotateFromTo(new double[] {0, 0, 1}, branchPos).getMatrix();
			for (int i = 0; i < numextrabranch; i++) {
				double[] extraPos;
				if (fibonacci) {
					double phi = acos(1 - 2 * (i + 0.5) / numextrabranch);
					double theta = PI * (1 + sqrt(5)) * i;
					double[] p = new double[] {cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi / 2), 1.0};
					Complex pos = stereographic(p);
					extraPos = inverseStereographic(pos);
				} else {
					double[] p = new double[] {rnd.nextGaussian(), rnd.nextGaussian(), Math.abs(rnd.nextGaussian()), 1.0};
					Pn.setToLength(p, p, 1, Pn.EUCLIDEAN);
					Complex pos = stereographic(p);
					pos = pos.times(pos);
					extraPos = inverseStereographic(pos);
				}
				T.transformVector(extraPos);
				CoVertex extraV = hds.addNewVertex();
				a.set(Position.class, extraV, extraPos);
			}
		}
		
		// equalization
		if (numEqualizerIterations > 0) {
			SphereUtility.equalizeSphereVertices(hds, branchVertices, numEqualizerIterations, 1E-6);
		}
		Set<CoEdge> glueSet = new HashSet<CoEdge>();
		int[] branchIndices = new int[branchVertices.size()];
		int i = 0;
		for (CoVertex bv : branchVertices) {
			branchIndices[i++] = bv.getIndex();
		}
		branchVertices = HyperellipticUtility.generateHyperellipticImage(hds, true, 0, true, glueSet, branchIndices);
		if (branchVerticesOUT != null) {
			branchVerticesOUT.addAll(branchVertices);
		}
		return hds;
	}
	
	@Override
	public void curveChanged(CurveChangeEvent e) {
		if (e.type == EventType.CURVE_CHANGED) {
			updatePeriodMatrix();
		}
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
		conformalDataPlugin = c.getPlugin(ConformalDataPlugin.class);
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