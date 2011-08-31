package de.varylab.discreteconformal.functional;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.halfedgetools.functional.Energy;
import de.jtem.halfedgetools.functional.Functional;
import de.varylab.discreteconformal.functional.FunctionalAdapters.InitialEnergy;

public interface ConformalFunctional <
	V extends Vertex<V, E, F>,
	E extends Edge<V, E, F>,
	F extends Face<V, E, F>
	> extends Functional<V, E, F> {

	public boolean triangleEnergyAndAlphas(
		// input	
			final DomainValue u, 
			final F f,
		// output
			final Energy E,
			final InitialEnergy<F> initialEnergy
	);
	
	public double getLambda(double length);
	public double getLength(double lambda);
	
	
	public double getNewLength(E e, DomainValue u);
	
}
