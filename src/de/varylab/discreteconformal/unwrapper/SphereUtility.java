package de.varylab.discreteconformal.unwrapper;

import java.util.Set;

import de.varylab.discreteconformal.functional.ElectrostaticSphereFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanApplication.TaoGradient;
import de.varylab.discreteconformal.unwrapper.numerics.CEuclideanApplication.TaoU;
import de.varylab.discreteconformal.unwrapper.numerics.ConformalEnergy;
import de.varylab.jpetsc.InsertMode;
import de.varylab.jpetsc.Vec;
import de.varylab.jtao.Tao;
import de.varylab.jtao.Tao.GetSolutionStatusResult;
import de.varylab.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.varylab.jtao.TaoApplication;

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
			TaoU u = new TaoU(x);
			TaoGradient G = new TaoGradient(g);
			ConformalEnergy E = new ConformalEnergy();
			fun.evaluate(hds, u, E, G, null);
			g.assemble();
			return E.get();
		}
		
	}

	
	public static void equalizeSphereVertices(CoHDS hds, Set<CoVertex> fixed, int steps, double tol) {
		Tao.Initialize();
		ElectrostaticSphereFunctional<CoVertex, CoEdge, CoFace> fun = new ElectrostaticSphereFunctional<CoVertex, CoEdge, CoFace>();
		fun.setFixedVertices(fixed);
		SphereOptimizationApplication app = new SphereOptimizationApplication(hds, fun);
		
		int n = fun.getDimension(hds);
		Vec u = new Vec(n);
		for (CoVertex v : hds.getVertices()) {
			u.setValue(v.getIndex() * 3 + 0, v.getPosition().get(0), InsertMode.INSERT_VALUES);
			u.setValue(v.getIndex() * 3 + 1, v.getPosition().get(1), InsertMode.INSERT_VALUES);
			u.setValue(v.getIndex() * 3 + 2, v.getPosition().get(2), InsertMode.INSERT_VALUES);
		}
		app.setInitialSolutionVec(u);
		
		Tao optimizer = new Tao(Tao.Method.LMVM);
		optimizer.setApplication(app);
		optimizer.setTolerances(0, 0, 0, 0);
		optimizer.setGradientTolerances(tol, tol, tol);
		optimizer.setMaximumIterates(steps);
		optimizer.solve();
		GetSolutionStatusResult stat = optimizer.getSolutionStatus();
		String status = stat.toString().replace("getSolutionStatus : ", "");
		System.out.println("optimization status ------------------------------------");
		System.out.println(status);
		for (CoVertex v : hds.getVertices()) {
			int i = v.getIndex() * 3;
			v.getPosition().set(0, u.getValue(i + 0));
			v.getPosition().set(1, u.getValue(i + 1));
			v.getPosition().set(2, u.getValue(i + 2));
			v.getPosition().normalize();
		}
	}
	
}
