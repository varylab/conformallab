package de.varylab.discreteconformal.logging;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.Logger;

import de.varylab.discreteconformal.ConformalLab;

public class LoggingUtility {

	public static void initLogging() {
		LogManager lm = LogManager.getLogManager();
		File localConf = new File("logging.properties");
		File localCustomConf = new File("logcustom.properties");
		InputStream confIn = null;
		InputStream propsIn = null;
		if (localConf.exists() || localCustomConf.exists()) {
			File confFile = null;
			if (localCustomConf.exists()) {
				confFile = localCustomConf;
			} else {
				confFile = localConf;
			}
			try {
				confIn = new FileInputStream(confFile);
				propsIn = new FileInputStream(confFile);
			} catch (FileNotFoundException e) {}
		} else {
			confIn = ConformalLab.class.getResourceAsStream("logging.properties");
			propsIn = ConformalLab.class.getResourceAsStream("logging.properties");			
		}
		assert confIn != null;
		assert propsIn != null;
		try {
			lm.reset();
			lm.readConfiguration(confIn);
			
			// create loggers from config file
			Properties p = new Properties();
			p.load(propsIn);
			for (Object o : p.keySet()) {
				String key = (String)o;
				if (key.endsWith(".level")) {
					String loggerName = key.substring(0, key.length() - 6);
					if (loggerName.isEmpty()) continue;
					Logger.getLogger(loggerName);
				}
			}
		} catch (Exception e) {
			System.out.println(e);
		} finally {
			try {
				confIn.close();
				propsIn.close();
			} catch (Exception e) {
				Logger l = lm.getLogger(LoggingUtility.class.getName());
				l.log(Level.WARNING, "could not close property streams", e);
			};
		}
	}

}
