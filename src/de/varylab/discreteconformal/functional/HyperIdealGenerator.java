package de.varylab.discreteconformal.functional;

import static java.lang.Math.PI;

import java.io.InputStream;
import java.util.List;

import javax.xml.bind.JAXBException;

import org.junit.Assert;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.TypedAdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.algorithm.subdivision.StellarLinear;
import de.jtem.halfedgetools.algorithm.triangulation.Triangulator;
import de.varylab.conformallab.data.DataIO;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.HalfedgeEmbedding;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class HyperIdealGenerator {

	public static CoHDS createLawsonSquareTiledWithBranchPoints() {
		TypedAdapterSet<double[]> a = new ConformalAdapterSet().querySet(double[].class);
		StellarLinear sl = new StellarLinear();
		CoHDS hds = new CoHDS();
		sl.execute(createLawsonSquareTiledBase(), hds, a);
		Assert.assertTrue(HalfEdgeUtils.isValidSurface(hds, true));
		Assert.assertEquals(2, HalfEdgeUtils.getGenus(hds));
		int index = 0;
		for (CoVertex v : hds.getVertices()) {
			if (index < 4) {
				v.setSolverIndex(index++);
			} else {
				v.setSolverIndex(-1);
			}
			v.setTheta(2 * PI);
		}
		for (CoEdge e : hds.getPositiveEdges()) {
			e.setSolverIndex(index);
			e.getOppositeEdge().setSolverIndex(index++);
			if (index < 17) {
				e.setTheta(PI);
				e.getOppositeEdge().setTheta(PI);
			} else {
				e.setTheta(PI/2);
				e.getOppositeEdge().setTheta(PI/2);
			}
		}
		return hds;
	}
	
	public static CoHDS createLawsonSquareTiled() {
		CoHDS hds = createLawsonSquareTiledBase();
		List<CoEdge> auxEdges = Triangulator.triangulateSingleSource(hds);
		Assert.assertTrue(HalfEdgeUtils.isValidSurface(hds, true));
		Assert.assertEquals(2, HalfEdgeUtils.getGenus(hds));
		int index = 0;
		for (CoVertex v : hds.getVertices()) {
			v.setSolverIndex(index++);
			v.setTheta(2 * PI);
		}
		for (CoEdge e : hds.getPositiveEdges()) {
			e.setSolverIndex(index);
			e.getOppositeEdge().setSolverIndex(index++);
			if (auxEdges.contains(e)) {
				e.setTheta(PI);
				e.getOppositeEdge().setTheta(PI);
			} else {
				e.setTheta(PI/2);
				e.getOppositeEdge().setTheta(PI/2);
			}
		}
		return hds;
	}

	
	private static CoHDS createLawsonSquareTiledBase() {
		CoHDS hds = new CoHDS(); 
		CoVertex A = hds.addNewVertex();
		CoVertex B = hds.addNewVertex();
		CoVertex C = hds.addNewVertex();
		CoVertex D = hds.addNewVertex();
		hds.addNewEdges(24);
		hds.getEdge(0).linkOppositeEdge(hds.getEdge(6));
		hds.getEdge(1).linkOppositeEdge(hds.getEdge(21));
		hds.getEdge(2).linkOppositeEdge(hds.getEdge(4));
		hds.getEdge(3).linkOppositeEdge(hds.getEdge(23));
		hds.getEdge(5).linkOppositeEdge(hds.getEdge(11));
		hds.getEdge(7).linkOppositeEdge(hds.getEdge(9));
		hds.getEdge(8).linkOppositeEdge(hds.getEdge(14));
		hds.getEdge(10).linkOppositeEdge(hds.getEdge(12));
		hds.getEdge(13).linkOppositeEdge(hds.getEdge(19));
		hds.getEdge(15).linkOppositeEdge(hds.getEdge(17));
		hds.getEdge(16).linkOppositeEdge(hds.getEdge(22));
		hds.getEdge(18).linkOppositeEdge(hds.getEdge(20));
		hds.getEdge(0).setTargetVertex(A);
		hds.getEdge(1).setTargetVertex(B);
		hds.getEdge(2).setTargetVertex(D);
		hds.getEdge(3).setTargetVertex(C);
		hds.getEdge(4).setTargetVertex(B);
		hds.getEdge(5).setTargetVertex(A);
		hds.getEdge(6).setTargetVertex(C);
		hds.getEdge(7).setTargetVertex(D);
		hds.getEdge(8).setTargetVertex(D);
		hds.getEdge(9).setTargetVertex(C);
		hds.getEdge(10).setTargetVertex(A);
		hds.getEdge(11).setTargetVertex(B);
		hds.getEdge(12).setTargetVertex(C);
		hds.getEdge(13).setTargetVertex(D);
		hds.getEdge(14).setTargetVertex(B);
		hds.getEdge(15).setTargetVertex(A);
		hds.getEdge(16).setTargetVertex(A);
		hds.getEdge(17).setTargetVertex(B);
		hds.getEdge(18).setTargetVertex(D);
		hds.getEdge(19).setTargetVertex(C);
		hds.getEdge(20).setTargetVertex(B);
		hds.getEdge(21).setTargetVertex(A);
		hds.getEdge(22).setTargetVertex(C);
		hds.getEdge(23).setTargetVertex(D);
		
		hds.getEdge(0).linkNextEdge(hds.getEdge(1));
		hds.getEdge(1).linkNextEdge(hds.getEdge(2));
		hds.getEdge(2).linkNextEdge(hds.getEdge(3));
		hds.getEdge(3).linkNextEdge(hds.getEdge(0));
		
		hds.getEdge(4).linkNextEdge(hds.getEdge(5));
		hds.getEdge(5).linkNextEdge(hds.getEdge(6));
		hds.getEdge(6).linkNextEdge(hds.getEdge(7));
		hds.getEdge(7).linkNextEdge(hds.getEdge(4));
		
		hds.getEdge(8).linkNextEdge(hds.getEdge(9));
		hds.getEdge(9).linkNextEdge(hds.getEdge(10));
		hds.getEdge(10).linkNextEdge(hds.getEdge(11));
		hds.getEdge(11).linkNextEdge(hds.getEdge(8));
		
		hds.getEdge(12).linkNextEdge(hds.getEdge(13));
		hds.getEdge(13).linkNextEdge(hds.getEdge(14));
		hds.getEdge(14).linkNextEdge(hds.getEdge(15));
		hds.getEdge(15).linkNextEdge(hds.getEdge(12));
		
		hds.getEdge(16).linkNextEdge(hds.getEdge(17));
		hds.getEdge(17).linkNextEdge(hds.getEdge(18));
		hds.getEdge(18).linkNextEdge(hds.getEdge(19));
		hds.getEdge(19).linkNextEdge(hds.getEdge(16));
		
		hds.getEdge(20).linkNextEdge(hds.getEdge(21));
		hds.getEdge(21).linkNextEdge(hds.getEdge(22));
		hds.getEdge(22).linkNextEdge(hds.getEdge(23));
		hds.getEdge(23).linkNextEdge(hds.getEdge(20));
		
		HalfEdgeUtils.fillAllHoles(hds);
		return hds;
	}

	public static CoHDS createLawsonHyperelliptic() throws JAXBException {
		InputStream in = HyperIdealGenerator.class.getResourceAsStream("lawson_curve_source.xml");
		HalfedgeEmbedding he = DataIO.readConformalData(HalfedgeEmbedding.class, in);
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<>();
		CoHDS hds = new CoHDS();
		AdapterSet a = new ConformalAdapterSet();
		DataUtility.toHalfedge(he, a, Position.class, hds, cutInfo);
		HyperIdealHyperellipticUtility.calculateCircleIntersections(hds);
		
		// vertex data
		int index = 0;
		for (CoVertex v : hds.getVertices()) {
			v.setTheta(2*PI);
			switch (v.getIndex()) {
			case 0: case 1: case 2: case 3: case 6: case 7:
				v.setSolverIndex(index++);
				break;
			default:
				v.setSolverIndex(-1);
				break;
			}
		}
		// edge angles and indices
		for (CoEdge e : hds.getPositiveEdges()) {
			e.setSolverIndex(index);
			e.getOppositeEdge().setSolverIndex(index++);
		}
		return hds;
	}

}
