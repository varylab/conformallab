package de.varylab.discreteconformal.logging;

import java.util.logging.Formatter;
import java.util.logging.LogRecord;

public class TrivialLogFormatter extends Formatter {

	@Override
	public String format(LogRecord log) {
		return log.getLevel() + ": " + log.getMessage();
	}

}
