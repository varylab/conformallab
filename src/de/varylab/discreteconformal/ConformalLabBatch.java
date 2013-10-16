package de.varylab.discreteconformal;

import joptsimple.OptionParser;
import joptsimple.OptionSet;

public class ConformalLabBatch {

	public void process(String... args) throws Exception {
		OptionParser p = new OptionParser();
		p.accepts("help", "Prints Help Information");
		
		OptionSet opts = null;
		try {
			opts = p.parse(args);
		} catch (Exception e) {
			System.out.println("Error: " + e.getMessage());
			printUsage(p);
			return;
		}
		if (opts.has("help")) {
			printUsage(p);
			return;
		}
		if (opts.nonOptionArguments().isEmpty()) {
			printUsage(p);
			return;
		}
	}
	
	
	private void printUsage(OptionParser p) throws Exception {
		System.out.println("Usage: java -jar DiscreteConformalLab.jar [Options] <data file>");
		p.printHelpOn(System.out);
	}
	
}
