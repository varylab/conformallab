package de.varylab.discreteconformal.plugin.hyperelliptic;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.util.List;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import de.jtem.java2d.DragListener;
import de.jtem.java2d.SceneComponent;
import de.jtem.java2d.TransformedMouseEvent;
import de.jtem.mfc.field.Complex;
import de.jtem.riemann.surface.BranchPoint;
import de.varylab.discreteconformal.plugin.hyperelliptic.CurveChangeEvent.EventType;

@SuppressWarnings("serial")
public class CurveEditor extends Editor implements CurveChangeListener {

	private Curve curve;
	private CoordinateSpinnerFrame coordFrame;

	private List<Point2D.Double> points;

	public CurveEditor(Curve curve) {
		super();

		initViewer();
		initCoordFrame();
		initMenu();

		setCurve(curve);

		update();
	}

	private void initCoordFrame() {
		coordFrame = new CoordinateSpinnerFrame();
		coordFrame.setSize(new Dimension(300, 40));
		coordFrame.setVisible(false);
	}

	private void initViewer() {
		addMouseListener(new MouseListener() {
			public void mouseClicked(MouseEvent arg0) {
				if (arg0.getButton() == MouseEvent.BUTTON2
						&& arg0.getClickCount() == 2) {
					double dx = arg0.getX(), dy = arg0.getY();
					updateUI();
					Rectangle2D.Double r = getViewport();
					double h = arg0.getComponent().getHeight();
					double w = arg0.getComponent().getWidth();
					dx /= w;
					dy /= h;
					dy = 1 - dy;
					double x = r.x + r.width * dx;
					double y = r.y + r.height * dy;
					addBrachPoint(x, y);
				}
			}

			public void mouseEntered(MouseEvent arg0) {
			}

			public void mouseExited(MouseEvent arg0) {
			}

			public void mousePressed(MouseEvent arg0) {
			}

			public void mouseReleased(MouseEvent arg0) {
			}
		});
	}

	private void removeBrachPoint(int id) {
		BranchPoint[] oldPoints = curve.getBranchPoints();
		BranchPoint[] newPoints = new BranchPoint[oldPoints.length - 1];
		int j = 0;
		for (int i = 0; i < oldPoints.length; i++) {
			if (i==id) {
				j--;
			} else {
				newPoints[i + j] = new BranchPoint(oldPoints[i].getCoords());
			}
		}
		setBranchPoints(newPoints);
	}

	private void addBrachPoint(double x, double y) {

		BranchPoint[] bPoints = new BranchPoint[curve.getNumOfBranchPoints() + 1];
		bPoints[0] = new BranchPoint(new Complex(x, y));
		for (int i = 1; i < bPoints.length; i++) {
			bPoints[i] = new BranchPoint(curve.getBranchPoint(i - 1)
					.getCoords());
		}

		setBranchPoints(bPoints);

	}

	private void setBranchPoints(BranchPoint[] bPoints) {
		Curve newCurve = new Curve(bPoints);
		newCurve.update();
		setCurve(newCurve);
	}

	protected void initScene() {

		scene = new SceneComponent();

		scene.setOutlinePaint(Color.BLACK);
		scene.setPointPaint(Color.RED);
		scene.setPointStroke(new BasicStroke(3));
		scene.setPointDragEnabled(true);
		scene.setPointOutlined(true);
		scene.setPointFilled(true);

		points = scene.getPoints();

		scene.addDragListener(new DragListener() {

			private double oldX, oldY;
			private int id, IDofCoordFrame;

			public void dragStart(TransformedMouseEvent e) {
				id = e.getIndexOfHitPart();
				double x = points.get(id).x;
				double y = points.get(id).y;
				if (e.getMouseEvent().getClickCount() != 2) {
					if (e.getMouseEvent().getButton() == MouseEvent.BUTTON2) {
						if (curve.getGenus() > 1)
							removeBrachPoint(id);
					} else {
						oldX = e.getX() - x;
						oldY = e.getY() - y;
					}
				} else {

					coordFrame.removeChangeListener();
					coordFrame.setTitle("Set Coordinates Of Point #" + id);
					IDofCoordFrame = id;

					ChangeListener listener = new ChangeListener() {
						final int ID = IDofCoordFrame;

						public void stateChanged(ChangeEvent e) {
							Complex z = new Complex((Double) coordFrame
									.getX_Spinner().getValue(),
									(Double) coordFrame.getY_Spinner()
											.getValue());

							curve.getBranchPoint(ID).setCoords(z);
							curve.update();
							update();
							getRoot().fireAppearanceChange();
						}
					};

					coordFrame.setCoords(curve.getBranchPoint(id).getCoords());

					coordFrame.setChangeListener(listener);

					coordFrame.pack();
					coordFrame.setVisible(true);

				}
			}

			public void drag(TransformedMouseEvent e) {
				if (e.getMouseEvent().getClickCount() == 2)
					return;
				if (id < 0) {
					System.err.println("branchpoint could not be found!");
					return;
				}
				Complex z = new Complex(e.getX() - oldX, e.getY() - oldY);

				curve.getBranchPoint(id).setCoords(z);

				if (coordFrame.isVisible() && id == IDofCoordFrame)
					coordFrame.setCoords(z);
				curve.update();
				update();
				getRoot().fireAppearanceChange();
				fireEditChangeEvent();
			}

			public void dragEnd(TransformedMouseEvent e) {
			}

		});

	}

	private void updatePoints() {

		points.clear();

		for (int j = 0; j < curve.getNumOfBranchPoints(); j++) {
			points.add(new Point2D.Double(curve.getBranchPoint(j).getRe(),
					curve.getBranchPoint(j).getIm()));
		}
	}

	public void update() {
		updatePoints();
		repaint();
	}

	public Curve getCurve() {
		return curve;
	}

	public void setCurve(Curve curve) {
		curve.update();
		boolean genusChanged = false;
		if (curve.equals(this.curve))
			return;
		if (this.curve == null) {
			curve.addCurveChangeListener(this);
			genusChanged = true;
			this.curve = curve;
		} else {
			curve.setCurveChangeListeners(this.curve.getCurveChangeListeners());
			if (this.curve.getGenus() != curve.getGenus())
				genusChanged = true;
			this.curve = curve;
		}
		update();
		if (genusChanged) {
			CurveChangeEvent e = new CurveChangeEvent(this.curve, this,
					EventType.GENUS_CHANGED);
			this.curve.fireCurveChangeEvent(e);
		}
		CurveChangeEvent e = new CurveChangeEvent(this.curve, this,
				EventType.CURVE_CHANGED);
		this.curve.fireCurveChangeEvent(e);
	}

	private JFileChooser fileChooser;

	private void initMenu() {

		initFileChooser();

		JMenuItem loadMenuItem = new JMenuItem("Load Surface");
		loadMenuItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O,
				ActionEvent.CTRL_MASK + ActionEvent.ALT_MASK));
		loadMenuItem.getAccessibleContext().setAccessibleDescription(
				"load surface");
		loadMenuItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFrame f = new JFrame();
				f.setAlwaysOnTop(true);
				f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				if (fileChooser.showOpenDialog(f) == JFileChooser.APPROVE_OPTION) {
					loadSurface(fileChooser.getSelectedFile());
				}
			}
		});

		JMenuItem saveMenuItem = new JMenuItem("Save Surface");
		saveMenuItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S,
				ActionEvent.CTRL_MASK + ActionEvent.ALT_MASK));
		saveMenuItem.getAccessibleContext().setAccessibleDescription(
				"save surface");
		saveMenuItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFrame f = new JFrame();
				f.setAlwaysOnTop(true);
				f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				if (fileChooser.showSaveDialog(f) == JFileChooser.APPROVE_OPTION) {
					saveSurface(fileChooser.getSelectedFile());
				}
			}
		});

		getMenu().add(loadMenuItem);
		getMenu().add(saveMenuItem);

	}

	private void initFileChooser() {
		fileChooser = new JFileChooser();
		fileChooser.addChoosableFileFilter(new CurveFileFilter());
		fileChooser.setAcceptAllFileFilterUsed(false);
		fileChooser.setFileView(new CurveFileView());
	}

	private void loadSurface(File f) {
		Curve t = CurveFileIO.readSpectralCurveFile(f);
		setCurve(t);
		CurveChangeEvent e = new CurveChangeEvent(this.curve, this,
				EventType.NEW_CURVE_SET);
		this.curve.fireCurveChangeEvent(e);
	}

	private void saveSurface(File file) {
		CurveFileIO.writeCurveToFile(getCurve(), file);
	}

	@Override
	public void curveChanged(CurveChangeEvent e) {
		if (e.type == EventType.CURVE_CHANGED)
			update();
	}

}
