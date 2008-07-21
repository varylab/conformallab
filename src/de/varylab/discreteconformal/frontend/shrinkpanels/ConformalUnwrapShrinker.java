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
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.math.optimization.NotConvergentException;
import de.varylab.discreteconformal.math.optimization.newton.NewtonOptimizer;
import de.varylab.discreteconformal.math.optimization.newton.NewtonOptimizer.Solver;

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


	public void widgetDefaultSelected(SelectionEvent e) {
		
	}


	public void widgetSelected(SelectionEvent e) {
		Object s = e.getSource();
		if (computeEnergyBtn == s) {
			computeEnergy();
		}
	}
	
	
	private void computeEnergy() {
		CHDS hds = ConformalLab.getGeometryController().getCHDS();
		DenseVector theta = new DenseVector(hds.numVertices());
		for (int i = 0; i < theta.size(); i++) {
			if (HalfEdgeUtils.isBoundaryVertex(hds.getVertex(i)))
				theta.set(i, 0.0);
			else
				theta.set(i, 2 * Math.PI);
		}
		hds.prepareInvariantData(theta);
		int n = hds.getDomainDimension();

		DenseVector u = new DenseVector(n);
		NewtonOptimizer optimizer = new NewtonOptimizer();
		optimizer.setSolver(Solver.GMRES);
		try {
			optimizer.minimize(u, hds);
		} catch (NotConvergentException e) {
			e.printStackTrace();
		}

		double[] E = {0.0};
		DenseVector G = new DenseVector(n);
		Matrix H = new DenseMatrix(n, n);
		hds.conformalEnergy(u, E, G, H);
		System.err.println("Energy: " + E[0]);
		System.err.println(G.toString());
		for (int i = 0; i < hds.numVertices(); i++) {
			for (int j = 0; j < hds.numVertices(); j++) {
				System.err.print(H.get(i, j) + "\t\t");
			}
			System.err.print("\n");
		}
	}
	
	
	
}
