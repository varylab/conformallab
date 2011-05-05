package de.varylab.discreteconformal.plugin;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.text.JTextComponent;

import de.jreality.plugin.basic.View;
import de.jreality.ui.LayoutFactory;
import de.jtem.blas.ComplexMatrix;
import de.jtem.jrworkspace.plugin.PluginInfo;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.surface.BranchPoint;
import de.jtem.riemann.theta.SiegelReduction;
import de.varylab.discreteconformal.hyperelliptic.Curve;
import de.varylab.discreteconformal.hyperelliptic.CurveChangeEvent;
import de.varylab.discreteconformal.hyperelliptic.CurveChangeEvent.EventType;
import de.varylab.discreteconformal.hyperelliptic.CurveChangeListener;
import de.varylab.discreteconformal.hyperelliptic.CurveEditor;
import de.varylab.discreteconformal.util.SimpleMatrixPrintUtility;

public class HyperellipticCurvePlugin extends ShrinkPanelPlugin implements
		CurveChangeListener {

	private JTextComponent matrixfield = new JTextArea();

	private CurveEditor editor;
	private JScrollPane protectorPane;

	public CurveEditor getEditor() {
		return editor;
	}

	public HyperellipticCurvePlugin() {

		Curve c = new Curve(1);
		c.update();

		setCurve(c);

		editor.setPreferredSize(new Dimension(300, 300));
		
		protectorPane = new JScrollPane(editor);
		protectorPane.setMinimumSize(editor.getPreferredSize());
		
		initViewMatrixPanel();

		setInitialPosition(SHRINKER_RIGHT);
		GridBagConstraints constraint = LayoutFactory.createRightConstraint();
		shrinkPanel.setTitle("Hyperelliptic Curve Plugin");
		constraint.weighty = 1.0;
		shrinkPanel.add(protectorPane, constraint);
		constraint.weighty = 0.0;
		shrinkPanel.add(matrixfield, constraint);

		shrinkPanel.setShrinked(false);
	}

	public Curve getCurve() {
		return editor.getCurve();
	}

	public void setCurve(Curve curve) {
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
	public void curveChanged(CurveChangeEvent e) {
		if (e.type == EventType.CURVE_CHANGED)
			updatePeriodMatrix();
	}

	private void initViewMatrixPanel() {
		matrixfield.setEditable(false);
		matrixfield.setText(" "
				+ SimpleMatrixPrintUtility.toString(
						getNormalizedPeriodMatrix(), 4));
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
		matrixfield.setText("" + SimpleMatrixPrintUtility.toString(
				getNormalizedPeriodMatrix(), 4));
		matrixfield.repaint();
	}

	@Override
	public PluginInfo getPluginInfo() {
		PluginInfo info = super.getPluginInfo();
		info.name = "Curve";
		info.vendorName = "";
		info.isDynamic = true;
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