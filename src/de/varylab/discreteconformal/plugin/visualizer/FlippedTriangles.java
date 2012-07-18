package de.varylab.discreteconformal.plugin.visualizer;

import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.plugin.data.DataSourceProvider;
import de.jtem.jrworkspace.plugin.Plugin;
import de.varylab.discreteconformal.adapter.TriangleOrientationAdapter;

public class FlippedTriangles extends Plugin implements DataSourceProvider {

	private TriangleOrientationAdapter adapter = new TriangleOrientationAdapter();	
	
	@Override
	public AdapterSet getDataSources() {
		return new AdapterSet(adapter);
	}

}
