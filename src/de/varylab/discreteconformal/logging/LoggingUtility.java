package de.varylab.discreteconformal.logging;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.logging.LogManager;

import de.varylab.discreteconformal.ConformalLab;

public class LoggingUtility {

	public static void initLogging() {
		LogManager lm = LogManager.getLogManager();
		File localConf = new File("logging.properties");
		File localCustomConf = new File("logcustom.properties");
		@SuppressWarnings("unused")
		String fileName = null;
		InputStream confIn = null;
		if (localConf.exists() | localCustomConf.exists()) {
			File confFile = null;
			if (localCustomConf.exists()) {
				confFile = localCustomConf;
			} else {
				confFile = localConf;
			}
			fileName = confFile.getAbsolutePath();
			try {
				confIn = new FileInputStream(confFile);
			} catch (FileNotFoundException e) {}
		} else {
			fileName = ConformalLab.class.getResource("logging.properties").getFile();
			confIn = ConformalLab.class.getResourceAsStream("logging.properties");
		}
		assert confIn != null;
		try {
			lm.readConfiguration(confIn);
		} catch (Exception e) {
			System.out.println(e);
		}
	}

}
