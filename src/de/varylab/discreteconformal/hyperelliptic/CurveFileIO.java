package de.varylab.discreteconformal.hyperelliptic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import de.jtem.blas.IntegerMatrix;
import de.jtem.blas.IntegerVector;
import de.jtem.riemann.surface.BranchPoint;

public class CurveFileIO {

	public static void writeCurveToFile(Curve curve, File file) {
		try {// try, catch - if something goes wrong

			// care for the right suffix

			String suffix = ".hc";

			File sf;
			int n = file.getPath().lastIndexOf(".");
			if (n < 0) {
				sf = new File(file.getPath() + suffix);
			} else {
				String str = file.getPath().substring(n,
						file.getPath().length());
				if (str.compareTo(suffix) == 0) {
					sf = file;
				} else {
					sf = new File(file.getPath() + suffix);
				}
			}

			// writing the document to write to
			PrintWriter out = new PrintWriter(new FileWriter(sf));

			// print
			out.println(getData(curve));

			// empty memory and close
			out.flush();
			out.close();
		} catch (Exception e) {
			System.err.println("couldn't save curve (exception caught): "
					+ e.getMessage());
		}

	}

	private static String getData(Curve curve) {
		String s = "";

		s += "genus " + curve.getGenus() + "\n";

		s += "\n";

		for (int i = 0; i < curve.getNumOfBranchPoints(); i++) {
			s += "point " + curve.getBranchPoint(i).getRe() + " "
					+ curve.getBranchPoint(i).getIm() + "\n";
		}

		s += "\n";

		// save transform history

		if (curve.getTransform().getG() != null) {
			s += "\n";
			s += "G " + curve.getTransform().getG().numRows + " "
					+ curve.getTransform().getG().numCols + " ";
			for (int i = 0; i < curve.getTransform().getG().numRows; i++) {
				for (int j = 0; j < curve.getTransform().getG().numCols; j++) {
					s += curve.getTransform().getG().re[i][j] + " ";
				}
			}
			s += "\n";
		}

		if (curve.getTransform().getGl() != null) {
			s += "\n";
			s += "Gl " + curve.getTransform().getGl().numRows + " "
					+ curve.getTransform().getGl().numCols + " ";
			for (int i = 0; i < curve.getTransform().getGl().numRows; i++) {
				for (int j = 0; j < curve.getTransform().getGl().numCols; j++) {
					s += curve.getTransform().getGl().re[i][j] + " ";
				}
			}
			s += "\n";
		}

		if (curve.getTransform().getT() != null) {
			s += "\n";
			s += "T " + curve.getTransform().getT().numRows + " "
					+ curve.getTransform().getT().numCols + " ";
			for (int i = 0; i < curve.getTransform().getT().numRows; i++) {
				for (int j = 0; j < curve.getTransform().getT().numCols; j++) {
					s += curve.getTransform().getT().re[i][j] + " ";
				}
			}
			s += "\n";
		}

		if (curve.getTransform().getH() != null) {
			s += "\n";
			s += "H " + curve.getTransform().getH().numRows + " "
					+ curve.getTransform().getH().numCols + " ";
			for (int i = 0; i < curve.getTransform().getH().numRows; i++) {
				for (int j = 0; j < curve.getTransform().getH().numCols; j++) {
					s += curve.getTransform().getH().re[i][j] + " ";
				}
			}
			s += "\n";
		}

		if (curve.getTransform().getDistinguishedPoints() != null) {
			s += "\n";
			s += "distinguishedPoints ";
			for (int i = 0; i < curve.getTransform().getDistinguishedPoints()
					.size(); i++) {
				s += curve.getTransform().getDistinguishedPoints().re[i] + " ";
			}
			s += "\n";
		}

		if (curve.getTransform().getMonodromyAboutInfinity() != null) {
			s += "\n";
			s += "monodromyAboutInfinity ";
			for (int i = 0; i < curve.getTransform()
					.getMonodromyAboutInfinity().size(); i++) {
				s += curve.getTransform().getMonodromyAboutInfinity().re[i]
						+ " ";
			}
			s += "\n";
		}

		if (curve.getTransform().getEdgeStartPoint() != null) {
			s += "\n";
			s += "egdeStartPoint ";
			for (int i = 0; i < curve.getTransform().getEdgeStartPoint().size(); i++) {
				s += curve.getTransform().getEdgeStartPoint().re[i] + " ";
			}
			s += "\n";
		}

		if (curve.getTransform().getEdgeBranchPoint() != null) {
			s += "\n";
			s += "egdeBranchPoint ";
			for (int i = 0; i < curve.getTransform().getEdgeBranchPoint()
					.size(); i++) {
				s += curve.getTransform().getEdgeBranchPoint().re[i] + " ";
			}
			s += "\n";
		}

		if (curve.getTransform().getEdgeEndPoint() != null) {
			s += "\n";
			s += "egdeEndPoint ";
			for (int i = 0; i < curve.getTransform().getEdgeEndPoint().size(); i++) {
				s += curve.getTransform().getEdgeEndPoint().re[i] + " ";
			}
			s += "\n";
		}

		s += "\n";
		s += "debug " + curve.getTransform().getDebug() + "\n";

		return s;
	}

	public static Curve readSpectralCurveFile(File f) {
		try {// open file, prepare to read
			FileReader fr = new FileReader(f);
			BufferedReader in = new BufferedReader(fr);
			try {

				if (in.markSupported() && in.ready()) {
					in.mark(9000);
				}

				// read the first line
				String line = in.readLine();

				// get data count
				int g = 0, pointcount = 0;

				// while there are lines
				while (line != null) {
					// and there is information in the line
					if (line.length() > 0) {
						// split the line by spaces
						String[] pieces = line.split(" ");
						if (pieces[0].compareTo("genus") == 0) {
							g = new Integer(pieces[1]);
						} else if (pieces[0].compareTo("point") == 0) {
							pointcount++;
						}
					}
					line = in.readLine();
				}

				in.reset();
				line = in.readLine();

				BranchPoint[] bPoints = new BranchPoint[pointcount];
				IntegerMatrix G = null, Gl = null, T = null, H = null;
				IntegerVector distinguishedPoints = null, monodromyAboutInfinity = null, edgeStartPoint = null, edgeBranchPoint = null, edgeEndPoint = null;
				boolean debug = false;
				int i = 0;

				// while there are lines
				while (line != null) {
					// and there is information in the line
					if (line.length() > 0) {
						// split the line by spaces
						String[] pieces = line.split(" ");
						if (pieces[0].compareTo("point") == 0) {
							bPoints[i] = new BranchPoint(new Double(pieces[1]),
									new Double(pieces[2]));
							i++;
						} else if (pieces[0].compareTo("G") == 0) {

							int numRows = new Integer(pieces[1]);
							int numCols = new Integer(pieces[2]);
							int[][] I = new int[numRows][numCols];
							for (int j = 0; j < numRows; j++) {
								for (int k = 0; k < numCols; k++) {
									I[j][k] = new Integer(pieces[j * numCols
											+ k + 3]);
								}
							}
							G = new IntegerMatrix(I);

						} else if (pieces[0].compareTo("Gl") == 0) {

							int numRows = new Integer(pieces[1]);
							int numCols = new Integer(pieces[2]);
							int[][] I = new int[numRows][numCols];
							for (int j = 0; j < numRows; j++) {
								for (int k = 0; k < numCols; k++) {
									I[j][k] = new Integer(pieces[j * numCols
											+ k + 3]);
								}
							}
							Gl = new IntegerMatrix(I);

						} else if (pieces[0].compareTo("T") == 0) {

							int numRows = new Integer(pieces[1]);
							int numCols = new Integer(pieces[2]);
							int[][] I = new int[numRows][numCols];
							for (int j = 0; j < numRows; j++) {
								for (int k = 0; k < numCols; k++) {
									I[j][k] = new Integer(pieces[j * numCols
											+ k + 3]);
								}
							}
							T = new IntegerMatrix(I);

						} else if (pieces[0].compareTo("H") == 0) {

							int numRows = new Integer(pieces[1]);
							int numCols = new Integer(pieces[2]);
							int[][] I = new int[numRows][numCols];
							for (int j = 0; j < numRows; j++) {
								for (int k = 0; k < numCols; k++) {
									I[j][k] = new Integer(pieces[j * numCols
											+ k + 3]);
								}
							}
							H = new IntegerMatrix(I);

						} else if (pieces[0].compareTo("distinguishedPoints") == 0) {
							int[] I = new int[pieces.length - 1];
							for (int j = 0; j < I.length; j++) {
								I[j] = new Integer(pieces[j + 1]);
							}
							distinguishedPoints = new IntegerVector(I);
						} else if (pieces[0]
								.compareTo("monodromyAboutInfinity") == 0) {
							int[] I = new int[pieces.length - 1];
							for (int j = 0; j < I.length; j++) {
								I[j] = new Integer(pieces[j + 1]);
							}
							monodromyAboutInfinity = new IntegerVector(I);
						} else if (pieces[0].compareTo("edgeStartPoint") == 0) {
							int[] I = new int[pieces.length - 1];
							for (int j = 0; j < I.length; j++) {
								I[j] = new Integer(pieces[j + 1]);
							}
							edgeStartPoint = new IntegerVector(I);
						} else if (pieces[0].compareTo("edgeBranchPoint") == 0) {
							int[] I = new int[pieces.length - 1];
							for (int j = 0; j < I.length; j++) {
								I[j] = new Integer(pieces[j + 1]);
							}
							edgeBranchPoint = new IntegerVector(I);
						} else if (pieces[0].compareTo("edgeEndPoint") == 0) {
							int[] I = new int[pieces.length - 1];
							for (int j = 0; j < I.length; j++) {
								I[j] = new Integer(pieces[j + 1]);
							}
							edgeEndPoint = new IntegerVector(I);
						} else if (pieces[0].equals("debug")) {
							debug = new Boolean(pieces[1]);
						}
					}
					line = in.readLine();
				}

				Curve curve = new Curve(g);
				curve.update();

				for (int j = 0; j < bPoints.length; j++) {
					curve.getBranchPoint(j).deconstrain();
					curve.getBranchPoint(j).setCoords(bPoints[j].getCoords());
				}

				in.reset();
				line = in.readLine();

				if (G != null)
					curve.getTransform().getG().assign(G);
				if (Gl != null)
					curve.getTransform().getGl().assign(Gl);
				if (T != null)
					curve.getTransform().getT().assign(T);
				if (H != null)
					curve.getTransform().getH().assign(H);
				if (distinguishedPoints != null)
					curve.getTransform().getDistinguishedPoints()
							.assign(distinguishedPoints);
				if (monodromyAboutInfinity != null)
					curve.getTransform().getMonodromyAboutInfinity()
							.assign(monodromyAboutInfinity);
				if (edgeStartPoint != null)
					curve.getTransform().getEdgeStartPoint()
							.assign(edgeStartPoint);
				if (edgeBranchPoint != null)
					curve.getTransform().getEdgeBranchPoint()
							.assign(edgeBranchPoint);
				if (edgeEndPoint != null)
					curve.getTransform().getEdgeEndPoint().assign(edgeEndPoint);

				curve.getTransform().setDebug(debug);

				curve.outdate();
				curve.update();

				// close the file
				in.close();

				return curve;

			} catch (IOException e) {
				System.err.println("couldn't load the curve! exception: "
						+ e.getMessage());
			}
		} catch (FileNotFoundException e) {

		}

		return null;

	}

}
