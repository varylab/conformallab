package de.varylab.discreteconformal.plugin.hyperelliptic;

import java.io.File;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.filechooser.FileView;

public class CurveFileView extends FileView {
	ImageIcon icon = new ImageIcon(getClass().getResource(
			"wenteTorus16.jpg"));

	@Override
	public String getName(File f) {
		return null; // let the L&F FileView figure this out
	}

	@Override
	public String getDescription(File f) {
		return null; // let the L&F FileView figure this out
	}

	@Override
	public Boolean isTraversable(File f) {
		return null; // let the L&F FileView figure this out
	}

	@Override
	public String getTypeDescription(File f) {
		String extension = Utils.getExtension(f);

		if (extension != null) {
			if (extension.equals(Utils.hc)) {
				return "Hyperelliptic Curve";
			}
		}
		return null;
	}

	@Override
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
