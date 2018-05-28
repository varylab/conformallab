package de.varylab.discreteconformal.unwrapper;

import java.util.Set;

import de.jreality.math.Pn;
import de.jtem.jpetsc.InsertMode;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.jtem.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.jtem.jtao.TaoApplication;
import de.varylab.discreteconformal.functional.ElectrostaticSphereFunctional;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Position;
import de.varylab.discreteconformal.functional.FunctionalAdapters.Variable;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.SimpleEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.TaoDomain;
import de.varylab.discreteconformal.unwrapper.numerics.TaoGradient;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CPosition;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;

public class SphereUtility {

	private static class SphereOptimizationApplication extends TaoApplication implements TaoAppAddCombinedObjectiveAndGrad {

		private ElectrostaticSphereFunctional<CoVertex, CoEdge, CoFace> 
			fun = null;
		private CoHDS
			hds = null;
		
		public SphereOptimizationApplication(
			CoHDS hds,
			ElectrostaticSphereFunctional<CoVertex, CoEdge, CoFace> fun
		) {
			super();
			this.hds = hds;
			this.fun = fun;
		}

		@Override
		public double evaluateObjectiveAndGradient(Vec x, Vec g) {
			TaoDomain u = new TaoDomain(x);
			TaoGradient G = new TaoGradient(g);
			SimpleEnergy E = new SimpleEnergy();
			fun.evaluate(hds, u, E, G, null);
			g.assemble();
			return E.get();
		}
		
	}

	
	public static void equalizeSphereVertices(CoHDS hds, Set<CoVertex> fixed, int steps, double tol) {
		Tao.Initialize();
		Variable<CoVertex, CoEdge> var = new CVariable();
		Position<CoVertex> pos = new CPosition();
		ElectrostaticSphereFunctional<CoVertex, CoEdge, CoFace> fun = new ElectrostaticSphereFunctional<CoVertex, CoEdge, CoFace>(var, pos);
		SphereOptimizationApplication app = new SphereOptimizationApplication(hds, fun);
		
		int index = 0; 
		for (CoVertex v : hds.getVertices()) {
			if (fixed.contains(v)) {
				v.setSolverIndex(-1);
			} else {
				v.setSolverIndex(index++);
			}
		}
		
		int n = fun.getDimension(hds);
		Vec u = new Vec(n);
		for (CoVertex v : hds.getVertices()) {
			if (!fixed.contains(v)) {
				index = v.getSolverIndex();
				Pn.dehomogenize(v.P, v.P);
				u.setValue(index * 3 + 0, v.P[0], InsertMode.INSERT_VALUES);
				u.setValue(index * 3 + 1, v.P[1], InsertMode.INSERT_VALUES);
				u.setValue(index * 3 + 2, v.P[2], InsertMode.INSERT_VALUES);
			}
		}
		app.setInitialSolutionVec(u);
		
		Tao optimizer = new Tao(Tao.Method.LMVM);
		optimizer.setApplication(app);
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setGradientTolerances(tol, tol, tol);
		optimizer.setMaximumIterates(steps);
		optimizer.solve();
		for (CoVertex v : hds.getVertices()) {
			if (!fixed.contains(v)) {
				int i = v.getSolverIndex() * 3;
				v.P[0] = u.getValue(i + 0);
				v.P[1] = u.getValue(i + 1);
				v.P[2] = u.getValue(i + 2);
				v.P[3] = 1.0;
			}
			Pn.setToLength(v.P, v.P, 1.0, Pn.EUCLIDEAN);
		}
	}
	
	
}
