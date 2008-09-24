package de.varylab.discreteconformal.frontend.shrinkpanels;

import static de.varylab.discreteconformal.ConformalLab.errorMessage;
import static de.varylab.discreteconformal.ConformalLab.getGeometryController;
import static org.eclipse.jface.layout.GridDataFactory.fillDefaults;
import static org.eclipse.swt.SWT.BORDER;
import static org.eclipse.swt.SWT.CHECK;
import static org.eclipse.swt.SWT.NONE;
import static org.eclipse.swt.SWT.SHADOW_ETCHED_IN;
import static org.eclipse.swt.layout.GridData.BEGINNING;
import static org.eclipse.swt.layout.GridData.CENTER;

import java.lang.reflect.InvocationTargetException;

import javax.swing.SwingUtilities;

import org.eclipse.core.runtime.IProgressMonitor;
import org.eclipse.jface.dialogs.ProgressMonitorDialog;
import org.eclipse.jface.layout.GridDataFactory;
import org.eclipse.jface.operation.IRunnableWithProgress;
import org.eclipse.swt.SWT;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.events.SelectionListener;
import org.eclipse.swt.layout.GridLayout;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.widgets.Group;
import org.eclipse.swt.widgets.Label;
import org.eclipse.swt.widgets.Spinner;

import de.varylab.discreteconformal.ConformalLab;
import de.varylab.discreteconformal.frontend.widget.ShrinkPanel;
import de.varylab.discreteconformal.frontend.widget.ShrinkPanelContainer;
import de.varylab.discreteconformal.heds.CHDS;
import de.varylab.discreteconformal.heds.unwrap.CDisk;
import de.varylab.discreteconformal.heds.unwrap.CSphere;
import de.varylab.discreteconformal.heds.unwrap.CUnwrapper;
import de.varylab.discreteconformal.heds.unwrap.UnwrapException;

public class UnwrapShrinker extends ShrinkPanel implements SelectionListener{

	private Button
		computeEnergyBtn = null;
	private Group
		coneConfigGroup = null;
	private Spinner
		numConesSpinner = null;
	private Button
		quantizeChecker = null;
	
	private int
		numCones = 0; 
	private boolean
		quantizeCones = true;
	
	
	public UnwrapShrinker(ShrinkPanelContainer parent) {
		super(parent, "Unwrap");
		createLayout();
		computeEnergyBtn.addSelectionListener(this);
		numConesSpinner.addSelectionListener(this);
		quantizeChecker.addSelectionListener(this);
	}

	
	private void createLayout() {
		setLayout(new GridLayout(1, true));
		
		coneConfigGroup = new Group(this, SHADOW_ETCHED_IN);
		coneConfigGroup.setText("Cone Singularities");
		coneConfigGroup.setLayout(new GridLayout(2, true));
		GridDataFactory.fillDefaults().grab(true, false).applyTo(coneConfigGroup);
		
		Label numConesLabel = new Label(coneConfigGroup, NONE);
		numConesLabel.setText("Cones");
		fillDefaults().grab(true, false).align(BEGINNING, CENTER).applyTo(numConesLabel);
		numConesSpinner = new Spinner(coneConfigGroup, BORDER);
		numConesSpinner.setValues(numCones, 0, 100, 0, 1, 1);
		fillDefaults().grab(true, false).applyTo(numConesSpinner);
		
		quantizeChecker = new Button(coneConfigGroup, CHECK);
		quantizeChecker.setText("Quantize Cone Angles");
		quantizeChecker.setSelection(true);
		fillDefaults().span(2,1).grab(true, false).applyTo(quantizeChecker);
		
		computeEnergyBtn = new Button(this, SWT.PUSH);
		computeEnergyBtn.setText("Unwrap");
		GridDataFactory.fillDefaults().grab(true, false).applyTo(computeEnergyBtn);
	}


	public void widgetDefaultSelected(SelectionEvent e) {
		
	}


	public void widgetSelected(SelectionEvent e) {
		Object s = e.getSource();
		if (computeEnergyBtn == s) {
		    try {
		        IRunnableWithProgress op = new ComputationThread();
		        new ProgressMonitorDialog(ConformalLab.getMainShell()).run(true, true, op);
		     } catch (InvocationTargetException ite) {
		    	 ite.printStackTrace();
		     } catch (InterruptedException ie) {
		    	 ie.printStackTrace();
		     }
		}
		if (numConesSpinner == s) {
			numCones = numConesSpinner.getSelection();
		}
		if (quantizeChecker == s) {
			quantizeCones = quantizeChecker.getSelection();
		}
	}
	
	
	private class ComputationThread implements IRunnableWithProgress {

		public void run(IProgressMonitor mon) throws InvocationTargetException, InterruptedException {
			final CHDS hds = getGeometryController().getCHDS();
			hds.setTexCoordinatesValid(false);
			
			// topology
			int X = hds.numVertices() - hds.numEdges() / 2 + hds.numFaces();
			CUnwrapper unwrapper = null;
			switch (X) {
				case 1:
					unwrapper = new CDisk(numCones, quantizeCones);
					break;
				case 2:
					unwrapper = new CSphere(numCones, quantizeCones);
					break;
				default:
					errorMessage("Error", "Unsupported topology");
					mon.setCanceled(true);
					return;
			}
			
			// unwrap
			try {
				unwrapper.unwrap(hds, mon);
			} catch (UnwrapException ue) {
				errorMessage("Error", ue.getMessage());
				mon.setCanceled(true);
				return;
			}
			
			hds.setTexCoordinatesValid(true);			
			// set geometry
			SwingUtilities.invokeLater(new Runnable(){
				public void run() {
					getGeometryController().setGeometry(hds);
					ConformalLab.getUIController().encompass();
				}
			});
		}
		
	}
	
}
