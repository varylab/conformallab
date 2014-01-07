package de.varylab.discreteconformal.heds;

import de.varylab.discreteconformal.functional.node.ConformalFace;

public class CoFace extends ConformalFace<CoVertex, CoEdge, CoFace> {

    public double[]
    	P = {0,0,0,1},
    	T = {0,0,0,1};

    @Override
    public void copyData(CoFace f) {
    	super.copyData(f);
    	P = f.P.clone();
    	T = f.T.clone();
    }
    
}
