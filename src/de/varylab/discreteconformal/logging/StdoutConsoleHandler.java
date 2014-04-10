package de.varylab.discreteconformal.logging;

import java.io.PrintStream;
import java.util.logging.ConsoleHandler;
import java.util.logging.Level;
import java.util.logging.LogRecord;

public class StdoutConsoleHandler extends ConsoleHandler {

	@Override
	public void publish(LogRecord record) {
		PrintStream out = System.out;
		if (record.getLevel().equals(Level.SEVERE) ||
			record.getLevel().equals(Level.WARNING)
		) {
			out = System.err;
		}
		out.println(getFormatter().format(record));
	}
	
	@Override
	public void close() {
		System.out.flush();
	}

}
