package de.varylab.discreteconformal.plugin.riemannmap;

import java.io.InputStream;

import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.varylab.conformallab.data.DataIO;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.HalfedgeMap;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoHDS;

public class CircularRiemannMap {

	private static InputStream
		dataIn = CircularRiemannMap.class.getResourceAsStream("square_noears.xml");
	private static AdapterSet 
		a = new ConformalAdapterSet();
	private static CoHDS 
		hds = new CoHDS();
	
	public static void main(String[] args) throws Exception {
		HalfedgeMap map = DataIO.readConformalData(HalfedgeMap.class, dataIn);
		DataUtility.toHalfedge(map.getImage(), a, Position.class, hds, null);
		System.out.println(hds);
		
	}
	
}
