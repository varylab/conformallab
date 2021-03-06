package de.varylab.discreteconformal.plugin.hyperelliptic;

import java.io.File;

import javax.swing.filechooser.FileFilter;

public class CurveFileFilter extends FileFilter {

	@Override
	public boolean accept(File f) {
		if (f.isDirectory()) {
			return true;
		}

		String extension = Utils.getExtension(f);
		if (extension != null) {
			if (extension.equals(Utils.hc)) {
				return true;
			} else {
				return false;
			}
		}

		return false;
	}

	// The description of this filter
	@Override
	public String getDescription() {
		return "Just Hyperelliptic Curves";
	}
}
