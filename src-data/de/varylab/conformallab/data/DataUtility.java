package de.varylab.conformallab.data;

import java.lang.annotation.Annotation;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

import de.jreality.geometry.IndexedFaceSetFactory;
import de.jreality.scene.IndexedFaceSet;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.varylab.conformallab.data.types.Circle;
import de.varylab.conformallab.data.types.DiscreteEmbedding;
import de.varylab.conformallab.data.types.DiscreteMetric;
import de.varylab.conformallab.data.types.EmbeddedTriangle;
import de.varylab.conformallab.data.types.EmbeddedVertex;
import de.varylab.conformallab.data.types.MetricEdge;
import de.varylab.conformallab.data.types.MetricTriangle;
import de.varylab.conformallab.data.types.ObjectFactory;
import de.varylab.conformallab.data.types.SchottkyData;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.MappedEdgeLengthAdapter;
import de.varylab.discreteconformal.plugin.schottky.SchottkyGenerator;

public class DataUtility {

	public static IndexedFaceSet toIndexedFaceSet(DiscreteEmbedding de) {
		int numVertices = de.getVertices().size();
		int numFaces = de.getTriangles().size();
	
		double[][] vData = new double[numVertices][];
		for (EmbeddedVertex ev : de.getVertices()) {
			vData[ev.getIndex()] = new double[]{ev.getX(), ev.getY(), ev.getZ(), ev.getW()};
		}
		int faceIndex = 0;
		int[][] fData = new int[numFaces][];
		for (EmbeddedTriangle et : de.getTriangles()) {
			fData[faceIndex++] = new int[]{et.getVertex1(), et.getVertex2(), et.getVertex3()};
		}
		
		IndexedFaceSetFactory f = new IndexedFaceSetFactory();
		f.setVertexCount(numVertices);
		f.setFaceCount(numFaces);
		f.setVertexCoordinates(vData);
		f.setFaceIndices(fData);
		f.setGenerateEdgesFromFaces(true);
		f.setGenerateFaceNormals(true);
		f.setGenerateVertexNormals(true);
		f.update();
		IndexedFaceSet ifs = f.getIndexedFaceSet();
		return ifs;
	}

	public static DiscreteMetric toDiscreteMetric(String name, CoHDS surface, AdapterSet a) {
		DiscreteMetric dm = new DiscreteMetric();
		dm.setName(name);
		int edgeIndex = 0;
		Map<Integer, Integer> edgeIndexMap = new HashMap<Integer, Integer>();
		for (CoEdge e : surface.getPositiveEdges()) {
			MetricEdge me = new MetricEdge();
			me.setIndex(edgeIndex);
			me.setLength(a.get(Length.class, e, Double.class));
			dm.getEdges().add(me);
			edgeIndexMap.put(e.getIndex(), edgeIndex);
			edgeIndex++;
		}
		for (CoFace f : surface.getFaces()) {
			MetricTriangle mt = new MetricTriangle();
			CoEdge e1 = f.getBoundaryEdge();
			CoEdge e2 = e1.getNextEdge();
			CoEdge e3 = e2.getNextEdge();
			e1 = e1.isPositive() ? e1 : e1.getOppositeEdge(); 
			e2 = e2.isPositive() ? e2 : e2.getOppositeEdge(); 
			e3 = e3.isPositive() ? e3 : e3.getOppositeEdge();
			mt.setEdge1(edgeIndexMap.get(e1.getIndex()));
			mt.setEdge2(edgeIndexMap.get(e2.getIndex()));
			mt.setEdge3(edgeIndexMap.get(e3.getIndex()));
			dm.getTriangles().add(mt);
		}
		return dm;
	}

	public static DiscreteEmbedding toDiscreteEmbedding(String name, CoHDS surface, AdapterSet a, Class<? extends Annotation> type) {
		DiscreteEmbedding de = new DiscreteEmbedding();
		de.setName(name);
		for (CoVertex v : surface.getVertices()) {
			EmbeddedVertex ev = new EmbeddedVertex();
			double[] pos = a.getD(type, v);
			ev.setX(pos[0]);
			ev.setY(pos[1]);
			ev.setZ(pos[2]);
			ev.setW(pos[3]);
			ev.setIndex(v.getIndex());
			de.getVertices().add(ev);
		}
		for (CoFace f : surface.getFaces()) {
			List<CoEdge> b = HalfEdgeUtils.boundaryEdges(f);
			assert b.size() == 3 : "only triangulations are supported";
			CoVertex v0 = b.get(0).getStartVertex();
			CoVertex v1 = b.get(1).getStartVertex();
			CoVertex v2 = b.get(2).getStartVertex();
			EmbeddedTriangle et = new EmbeddedTriangle();
			et.setVertex1(v0.getIndex());
			et.setVertex2(v1.getIndex());
			et.setVertex3(v2.getIndex());
			de.getTriangles().add(et);
		}
		return de;
	}

	public static CoHDS toHalfedgeAndLengths(DiscreteMetric dm, MappedEdgeLengthAdapter lengthAdapterOutput) {
		CoHDS hds = new CoHDS();
		Map<MetricEdge, CoEdge> edgeMap = new HashMap<MetricEdge, CoEdge>();
		for (MetricEdge me : dm.getEdges()) {
			CoEdge e1 = hds.addNewEdge();
			CoEdge e2 = hds.addNewEdge();
			e1.linkOppositeEdge(e2);
			edgeMap.put(me, e1);
		}
		Set<CoEdge> usedEdges = new HashSet<CoEdge>();
		for (MetricTriangle mt : dm.getTriangles()) {
			CoFace f = hds.addNewFace();
			MetricEdge e1 = dm.getEdges().get(mt.getEdge1());
			MetricEdge e2 = dm.getEdges().get(mt.getEdge2());
			MetricEdge e3 = dm.getEdges().get(mt.getEdge3());
			CoEdge ce1 = edgeMap.get(e1);
			CoEdge ce2 = edgeMap.get(e2);
			CoEdge ce3 = edgeMap.get(e3);
			if (usedEdges.contains(ce1)) {
				ce1 = ce1.getOppositeEdge();
			}
			if (usedEdges.contains(ce2)) {
				ce2 = ce2.getOppositeEdge();
			}
			if (usedEdges.contains(ce3)) {
				ce3 = ce3.getOppositeEdge();
			}
			usedEdges.add(ce1);
			usedEdges.add(ce2);
			usedEdges.add(ce3);
			ce1.linkNextEdge(ce2);
			ce2.linkNextEdge(ce3);
			ce3.linkNextEdge(ce1);
			ce1.setLeftFace(f);
			ce2.setLeftFace(f);
			ce3.setLeftFace(f);
		}
		// link boundary
		Set<CoEdge> boundaryEdges = new HashSet<CoEdge>(hds.getEdges());
		boundaryEdges.removeAll(usedEdges);
		for (CoEdge be : boundaryEdges) {
			CoEdge e = be.getOppositeEdge();
			CoEdge ne = e.getPreviousEdge().getOppositeEdge();
			while (ne.getLeftFace() != null) {
				ne = ne.getPreviousEdge().getOppositeEdge();
			}
			be.linkNextEdge(ne);
		}
		Queue<CoEdge> linkQueue = new LinkedList<CoEdge>(hds.getEdges());
		while (!linkQueue.isEmpty()) {
			CoEdge e = linkQueue.poll();
			if (e.getTargetVertex() != null) {
				continue;
			}
			CoVertex v = hds.addNewVertex();
			CoEdge ne = e;
			do {
				ne.setTargetVertex(v);
				ne = ne.getNextEdge().getOppositeEdge();
			} while (ne != e);
		}
		
		for (MetricEdge e : dm.getEdges()) {
			CoEdge e1 = edgeMap.get(e);
			CoEdge e2 = e1.getOppositeEdge();
			lengthAdapterOutput.setEdgeValue(e1, e.getLength(), null);
			lengthAdapterOutput.setEdgeValue(e2, e.getLength(), null);
		}
		return hds;
	}

	public static SchottkyData toSchottkyData(String name, List<SchottkyGenerator> generators) {
		ObjectFactory of = new ObjectFactory();
		SchottkyData schottkyData = of.createSchottkyData();
		schottkyData.setName(name);
		for (SchottkyGenerator s : generators) {
			de.varylab.conformallab.data.types.SchottkyGenerator g = of.createSchottkyGenerator();
			de.varylab.conformallab.data.types.Complex Ad = of.createComplex();
			de.varylab.conformallab.data.types.Complex Bd = of.createComplex();
			de.varylab.conformallab.data.types.Complex Mud = of.createComplex();
			de.varylab.conformallab.data.types.Complex center = of.createComplex();
			Circle Cd = of.createCircle();
			Ad.setRe(s.getA().re);
			Ad.setIm(s.getA().im);
			Bd.setRe(s.getB().re);
			Bd.setIm(s.getB().im);
			Mud.setRe(s.getMu().re);
			Mud.setIm(s.getMu().im);
			center.setRe(s.getCycle().getCenter().re);
			center.setIm(s.getCycle().getCenter().im);
			Cd.setCenter(center);
			Cd.setRadius(s.getCycle().getRadius());
			g.setA(Ad);
			g.setB(Bd);
			g.setMu(Mud);
			g.setCircle(Cd);
			schottkyData.getGenerators().add(g);
		}
		return schottkyData;
	}

	public static List<SchottkyGenerator> toGeneratorsList(SchottkyData data) {
		List<SchottkyGenerator> result = new LinkedList<SchottkyGenerator>();
		for (de.varylab.conformallab.data.types.SchottkyGenerator s : data.getGenerators()) {
			SchottkyGenerator g = new SchottkyGenerator();
			g.getA().setRe(s.getA().getRe());
			g.getA().setIm(s.getA().getIm());
			g.getB().setRe(s.getB().getRe());
			g.getB().setIm(s.getB().getIm());
			g.getMu().setRe(s.getMu().getRe());
			g.getMu().setIm(s.getMu().getIm());
			g.getCycle().getCenter().setRe(s.getCircle().getCenter().getRe());
			g.getCycle().getCenter().setIm(s.getCircle().getCenter().getIm());
			g.getCycle().setRadius(s.getCircle().getRadius());
			result.add(g);
		}
		return result;
	}

}
