package de.varylab.discreteconformal.functional.node;

import de.jtem.halfedge.Face;

public class ConformalFace <
	V extends ConformalVertex<V, E, F>,
	E extends ConformalEdge<V, E, F>,
	F extends ConformalFace<V, E, F>
> extends Face<V, E, F> {

    private double 
    	initialEnergy = 0.0;
    
  	public double getInitialEnergy() {
		return initialEnergy;
	}

	public void setInitialEnergy(double energy) {
		this.initialEnergy = energy;
	}
}
