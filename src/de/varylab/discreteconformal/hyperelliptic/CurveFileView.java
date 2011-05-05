package de.varylab.discreteconformal.hyperelliptic;

import java.io.File;
import javax.swing.*;
import javax.swing.filechooser.*;

public class CurveFileView extends FileView {
	ImageIcon icon = new ImageIcon(getClass().getResource(
			"wenteTorus16.jpg"));

	public String getName(File f) {
		return null; // let the L&F FileView figure this out
	}

	public String getDescription(File f) {
		return null; // let the L&F FileView figure this out
	}

	public Boolean isTraversable(File f) {
		return null; // let the L&F FileView figure this out
	}

	public String getTypeDescription(File f) {
		String extension = Utils.getExtension(f);

		if (extension != null) {
			if (extension.equals(Utils.hc)) {
				return "Hyperelliptic Curve";
			}
		}
		return null;
	}

	public Icon getIcon(File f) {
		String extension = Utils.getExtension(f);

		if (extension != null) {
			if (extension.equals(Utils.hc)) {
				return icon;
			}
		}
		return null;
	}
}
