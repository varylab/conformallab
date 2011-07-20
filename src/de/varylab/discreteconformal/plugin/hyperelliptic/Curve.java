package de.varylab.discreteconformal.plugin.hyperelliptic;

import java.util.List;
import java.util.Vector;

import de.jtem.mfc.field.Complex;
import de.jtem.riemann.constrainedComplex.Fix;
import de.jtem.riemann.surface.BranchPoint;
import de.jtem.riemann.surface.DistinguishedPoint;
import de.jtem.riemann.surface.Origin;
import de.jtem.riemann.surface.hyperElliptic.HyperEllipticCurve;
import de.varylab.discreteconformal.plugin.hyperelliptic.CurveChangeEvent.EventType;

public class Curve extends HyperEllipticCurve {

	private static final long serialVersionUID = 1L;

	private List<CurveChangeListener> listeners = new Vector<CurveChangeListener>();
	
	public void addCurveChangeListener(
			CurveChangeListener listener) {
		listeners.add(listener);
	}

	public void removeCurveChangeListener(
			CurveChangeListener listener) {
		listeners.remove(listener);
	}

	public void setCurveChangeListeners(
			List<CurveChangeListener> listeners) {
		this.listeners = listeners;
	}

	public List<CurveChangeListener> getCurveChangeListeners() {
		return listeners;
	}

	public void fireCurveChangeEvent(CurveChangeEvent e) {
		for (CurveChangeListener l : getCurveChangeListeners()) {
			l.curveChanged(e);
		}
	}

	protected static BranchPoint[] getInitBranchPoints(int genus) {

		BranchPoint[] branchPoint = new BranchPoint[2 * genus + 2];

		branchPoint[0] = new BranchPoint(0, 1);
		branchPoint[1] = new BranchPoint(0, -1);

		for (int i = 0, I = 2; i < genus; i++, I += 2) {

			branchPoint[I] = new BranchPoint(i + 1., 1);
			branchPoint[I + 1] = new BranchPoint(i + 1., -1);

		}

		return branchPoint;
	}

	protected static DistinguishedPoint[] getInitDistinguishedPoints() {

		DistinguishedPoint[] distinguishedPoint = new DistinguishedPoint[1];

		new Fix(distinguishedPoint[0] = new DistinguishedPoint(0, 0));

		return distinguishedPoint;
	}

	static protected Complex getInitOrigin() {
		return new Complex(0.5, 0);
	}

	Curve() {
		distinguishedPoint = getInitDistinguishedPoints();
		origin = new Origin(0.5, 0);
	}

	public Curve(int aGenus) {

		super(getInitBranchPoints(aGenus), null, getInitDistinguishedPoints(),
				getInitOrigin());
	}

	public Curve(BranchPoint[] branchPoints) {
		super(branchPoints, null, getInitDistinguishedPoints(), getInitOrigin());
		setSymmetrizePeriodMatrix(false);
	}

	public void update() {
		super.update();
		CurveChangeEvent e = new CurveChangeEvent(this, this,
				EventType.CURVE_CHANGED);
		fireCurveChangeEvent(e);
	}

}
