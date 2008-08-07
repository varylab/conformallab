package de.varylab.discreteconformal.image;

import java.io.IOException;
import java.io.InputStream;

import javax.imageio.ImageIO;

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


	public static Image getSWTImage(String filename){
		InputStream in = ImageHook.class.getResourceAsStream(filename);
		if (in == null)
			return null;
		return  new Image(Display.getCurrent(), in);
	}
	
	
	public static java.awt.Image getAWTImage(String filename){
		InputStream in = ImageHook.class.getResourceAsStream(filename);
		if (in == null)
			return null;
		java.awt.Image result = null;
		try {
			result = ImageIO.read(in);
		} catch (IOException e) {}
		return result;
	}
	
}
