package de.varylab.discreteconformal.logging;

import java.util.logging.Formatter;
import java.util.logging.LogRecord;

public class SequenceNumberLogFormatter extends Formatter {

	@Override
	public String format(LogRecord log) {
		String className = log.getSourceClassName();
		className = className.substring(className.lastIndexOf('.') + 1);
		String methodName = log.getSourceMethodName();
		return "["+ log.getSequenceNumber() + "]" + className + "." +  methodName + "(): " + log.getMessage();
	}

}
