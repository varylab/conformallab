package de.varylab.discreteconformal.holomorphicformsexperiments;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;

/**
 * @author Kristoffer Josefsson, Andre Heydt 
 * TODO: Fix for meshes with boundary
 */
public class Loop {

	//return new HEDS approximated using dyadic scheme
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> Map<E, Set<E>> subdivide(
		HDS oldh, 
		HDS newh, 
		AdapterSet a
	){
		Map<V,E> newVtoOldE = new HashMap<V,E>();
		Map<E, Set<E>> oldEtoNewEs = new HashMap<E,Set<E>>();
		oldh.createCombinatoriallyEquivalentCopy(newh);

		
		for(E oe : oldh.getPositiveEdges()){
			Double l = a.get(Length.class, oe, Double.class);
			E e = newh.getEdge(oe.getIndex());

			V mv = newh.addNewVertex();
			newVtoOldE.put(mv, oe);
			
			V et = e.getTargetVertex();
			V es = e.getStartVertex();
			E eo = e.getOppositeEdge();
			
			E e2 = newh.addNewEdge();
			E eo2 = newh.addNewEdge();	
			E en = e.getNextEdge();
			E eon = eo.getNextEdge();		
			
		//	re-link edges	
			e.setTargetVertex(mv);
			e.linkNextEdge(e2);
			e.linkOppositeEdge(eo2);
			e2.linkNextEdge(en);
			e2.setTargetVertex(et);
			eo.setTargetVertex(mv);
			eo.linkNextEdge(eo2);
			eo.linkOppositeEdge(e2);
			eo2.linkNextEdge(eon);
			eo2.setTargetVertex(es);
			
			Set<E> newEs = new HashSet<E>();
			newEs.add(e); 
			newEs.add(e2);
			
			oldEtoNewEs.put(oe, newEs);

			a.set(Length.class, e, l / 2);
			a.set(Length.class, e.getOppositeEdge(), l / 2);
			a.set(Length.class, e2, l / 2);
			a.set(Length.class, eo2.getOppositeEdge(), l / 2);
		}
		
		
	//	end : edge cut
		
	//	rearrange interior
		for(F of : oldh.getFaces()){
			F f = newh.getFace(of.getIndex());
			List<E> e = new ArrayList<E>(0);
			List<E> eb = new ArrayList<E>(0);
			List<F> fn = newh.addNewFaces(3);
			
			e.add(f.getBoundaryEdge());

			eb.add(e.get(0).getNextEdge());
			e.add(eb.get(0).getNextEdge());
			eb.add(e.get(1).getNextEdge());
			e.add(eb.get(1).getNextEdge());
			eb.add(e.get(2).getNextEdge());
			
			List<E> inner = newh.addNewEdges(3);
			List<E> outer = newh.addNewEdges(3);
			
			for (int i=0 ; i<3; i++){
				inner.get(i).setLeftFace(f);
				inner.get(i).setTargetVertex(e.get(i).getTargetVertex());
				inner.get(i).linkOppositeEdge(outer.get(i));
				inner.get(i).linkNextEdge(inner.get((i+1)%3));
				
				e.get(i).linkNextEdge(outer.get(i));
				outer.get(i).linkNextEdge(eb.get((i+2)%3));
				outer.get(i).setTargetVertex(e.get((i+2)%3).getTargetVertex());
				e.get(i).setLeftFace(fn.get(i));
				outer.get(i).setLeftFace(fn.get(i));
				eb.get((i+2)%3).setLeftFace(fn.get(i));
				
			};
			
			for (E ee : inner) {
				E eee = ee.getPreviousEdge().getOppositeEdge().getNextEdge();
				Double l = a.get(Length.class, eee, Double.class);
				a.set(Length.class, ee, l);
				a.set(Length.class, ee.getOppositeEdge(), l);
			}
		}
		
		return oldEtoNewEs;
	};
	
}
