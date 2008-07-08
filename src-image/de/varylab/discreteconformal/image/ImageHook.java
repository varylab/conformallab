package de.varylab.discreteconformal.image;
import java.io.InputStream;

import org.eclipse.swt.graphics.Image;
import org.eclipse.swt.widgets.Display;

/**
 * An image loader for jar files
 * <p>
 * Copyright 2005 <a href="http://www.sechel.de">Stefan Sechelmann</a>
 * <a href="http://www.math.tu-berlin.de/geometrie">TU-Berlin</a> 
 * @author Stefan Sechelmann
 */
public class ImageHook {


	public static Image getImage(String filename){
		InputStream in = ImageHook.class.getResourceAsStream(filename);
		if (in == null)
			return null;
		return  new Image(Display.getCurrent(), in);
	}
	
	
}
