package de.varylab.discreteconformal.frontend.shrinkpanels;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;

import org.eclipse.jface.layout.GridDataFactory;
import org.eclipse.swt.SWT;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.events.SelectionListener;
import org.eclipse.swt.layout.GridLayout;
import org.eclipse.swt.widgets.Button;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.ConformalLab;
import de.varylab.discreteconformal.frontend.widget.ShrinkPanel;
import de.varylab.discreteconformal.frontend.widget.ShrinkPanelContainer;
import de.varylab.discreteconformal.heds.CEdge;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.CVertex;

public class ConformalUnwrapShrinker extends ShrinkPanel implements SelectionListener{

	private Button
		computeEnergyBtn = null;
	
	
	public ConformalUnwrapShrinker(ShrinkPanelContainer parent) {
		super(parent, "Conformal Unwrap");
		createLayout();
		computeEnergyBtn.addSelectionListener(this);
	}

	
	private void createLayout() {
		setLayout(new GridLayout(1, true));
		computeEnergyBtn = new Button(this, SWT.PUSH);
		computeEnergyBtn.setText("Compute Conformal Energy");
		GridDataFactory.fillDefaults().grab(true, false).applyTo(computeEnergyBtn);
	}


	@Override
	public void widgetDefaultSelected(SelectionEvent e) {
		
	}


	@Override
	public void widgetSelected(SelectionEvent e) {
		Object s = e.getSource();
		if (computeEnergyBtn == s) {
			computeEnergy();
		}
	}
	
	
	private void computeEnergy() {
		CHDS hds = ConformalLab.getGeometryController().getCHDS();
		DenseVector u = new DenseVector(hds.numVertices());
		DenseVector theta = new DenseVector(hds.numVertices());
		DenseVector G = new DenseVector(hds.numVertices());
		Matrix H = new DenseMatrix(hds.numVertices(), hds.numVertices());
		for (int i = 0; i < theta.size(); i++) {
			if (HalfEdgeUtils.isBoundaryVertex(hds.getVertex(i)))
				theta.set(i, 0.0);
			else
				theta.set(i, 2 * Math.PI);
		}
		double[] E = {0.0};
		hds.prepareInvariantData(theta);
		hds.conformalEnergy(u, E, G, H);
		System.err.println("Energy: " + E[0]);
		System.err.println(G.toString());
		for (int i = 0; i < hds.numVertices(); i++) {
			for (int j = 0; j < hds.numVertices(); j++) {
				System.err.print(H.get(i, j) + "\t\t");
			}
			System.err.print("\n");
		}
		for (CEdge e : hds.getPositiveEdges()) {
			System.err.println("Lambda: " + e.getLambda());
		}
		System.err.println("Vertices -------------------------");
		for (CVertex v: hds.getVertices()) 
			System.err.println("V" + v.getIndex() + ": "+ v.getPosition());
	}
	
	
	
}
