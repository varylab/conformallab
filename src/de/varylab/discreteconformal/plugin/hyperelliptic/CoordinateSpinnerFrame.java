package de.varylab.discreteconformal.plugin.hyperelliptic;

import java.awt.Dimension;
import java.awt.GridBagLayout;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeListener;

import de.jtem.mfc.field.Complex;

public class CoordinateSpinnerFrame extends JFrame{
	
	static final long serialVersionUID= 1L;
	private JSpinner X_Spinner, Y_Spinner;
	private ChangeListener listener= null;
	
	public CoordinateSpinnerFrame() {
		super("Set Coordinates");
		setAlwaysOnTop(true);
		setResizable(false);
		JPanel pane= new JPanel();
		pane.setPreferredSize(new Dimension(310,30));
		pane.setMaximumSize(new Dimension(310,30));
		pane.setMinimumSize(new Dimension(310,30));
		pane.setLayout(new GridBagLayout());
		
		pane.add(new JLabel("x: "));
		
		X_Spinner= new JSpinner(new SpinnerNumberModel(0.001,-1E+10,+1E+10,0.001));
		X_Spinner.setPreferredSize(new Dimension(130,30));
		X_Spinner.setMaximumSize(new Dimension(130,30));
		X_Spinner.setMinimumSize(new Dimension(130,30));
		pane.add(X_Spinner);
		
		pane.add(new JLabel("   y:"));
		
		Y_Spinner= new JSpinner(new SpinnerNumberModel(0.001,-1E+10,+1E+10,0.001));
		Y_Spinner.setPreferredSize(new Dimension(130,30));
		Y_Spinner.setMaximumSize(new Dimension(130,30));
		Y_Spinner.setMinimumSize(new Dimension(130,30));
		pane.add(Y_Spinner);
					
		add(pane);

	}
	
	public void setChangeListener(ChangeListener listener){
		this.listener= listener;
		X_Spinner.addChangeListener(listener);
		Y_Spinner.addChangeListener(listener);
	}

	public void removeChangeListener(){
		X_Spinner.removeChangeListener(listener);
		Y_Spinner.removeChangeListener(listener);
	}
	
	public JSpinner getX_Spinner() {
		return X_Spinner;
	}
	
	public void setX_Spinner(JSpinner spinner) {
		X_Spinner = spinner;
	}
	
	public JSpinner getY_Spinner() {
		return Y_Spinner;
	}
	
	public void setY_Spinner(JSpinner spinner) {
		Y_Spinner = spinner;
	}
	
	public void setCoords(Complex z){
		X_Spinner.setValue(z.re);
		Y_Spinner.setValue(z.im);
	}

	public ChangeListener getListener() {
		return listener;
	}
	
}
