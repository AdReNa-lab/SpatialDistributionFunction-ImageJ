//Spatial Distribution Function
//Copyright (C) 2019  Niamh Mac Fhionnlaoich, Runzhang Qi, Stefan Guldin
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <https://www.gnu.org/licenses/>.

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GUI;
import ij.gui.GenericDialog;
import ij.gui.MultiLineLabel;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.frame.PlugInFrame;
import ij.process.ImageProcessor;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.PolarChartPanel;
import org.jfree.chart.plot.PolarAxisLocation;
import org.jfree.chart.plot.PolarPlot;
import org.jfree.chart.renderer.DefaultPolarItemRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import utils.GUIHelper;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

public class SpatialDistributionFunction_Polar implements PlugInFilter {

	ImagePlus imp;
	private double RMinDefault;
	private double RMaxDefault;
	private double SpacingDefault;

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		return DOES_ALL+NO_UNDO;
	}
	@SuppressWarnings("Duplicates")
	public void run(ImageProcessor ip) {
		ResultsTable coordinateTable = Analyzer.getResultsTable();
		// Terminate here if table is not available. otherwise error will occur
		int xColumnIndex = coordinateTable.getColumnIndex("X");
		int yColumnIndex = coordinateTable.getColumnIndex("Y");
		// Read in positional information
		float[] xValues = coordinateTable.getColumn(xColumnIndex);
		float[] yValues = coordinateTable.getColumn(yColumnIndex);
		// Read from dialog
		double Rmin;
		double Rmax;
		float spacing;
		GenericDialog gd = new GenericDialog("2D RDF", IJ.getInstance());
		gd.addNumericField("R min (pixels):", RMaxDefault, 2);
		gd.addNumericField("R max (pixels):", RMinDefault, 2);
		gd.addNumericField("Spacing", SpacingDefault, 2);
		gd.addMessage("This plugin analyses colloidal ordering via the SDF and plots the resulting analysis in Cartesian and polar coordinates. \n" +
				"(see scientific publication in Langmuir:  https://doi.org/10.1021/acs.langmuir.9b02877)");
		MultiLineLabel text = (MultiLineLabel) gd.getMessage();
		GUIHelper.addHyperLinkListener(text, "https://doi.org/10.1021/acs.langmuir.9b02877");
		gd.showDialog();
		if (gd.wasCanceled()) {
			return;
		}
		RMinDefault = gd.getNextNumber();
		RMaxDefault = gd.getNextNumber();
		SpacingDefault = gd.getNextNumber();
		Rmin = (int) RMinDefault;
		Rmax = (int) RMaxDefault;
		spacing = (float) SpacingDefault;
		if (gd.invalidNumber()) {
			IJ.showMessage("Error", "Invalid input Number");
			return;
		}

		float[] Xd = xValues.clone();
		float[] Yd = yValues.clone();
		// Pre-sets the Xa, Ya vectors for speed
		float[] Xa = new float[Xd.length * Xd.length];
		float[] Ya = new float[Yd.length * Yd.length];
		// Create array for every particle in relation to every other particle
		System.out.println("2D RDF: Create array for every particle in relation to every other particle...");
		for (int i = 0; i < Xd.length; i++) {
			// Chooses the ith particle as the reference particle
			float Xc = Xd[i];
			float Yc = Yd[i];
			// Shifts the array such that the reference particle lies at the 0,0
			// point
			float[] X_temp = new float[Xd.length];
			float[] Y_temp = new float[Yd.length];
			for (int j = 0; j < Xd.length; j++) {
				X_temp[j] = Xd[j] - Xc;
				Y_temp[j] = Yd[j] - Yc;
			}
			// Adds the above array to the Xa, Ya vectors
			for (int k = 0; k < Xd.length; k++) {
				Xa[i * Xd.length + k] = X_temp[k];
				Ya[i * Xd.length + k] = Y_temp[k];
			}
		}

		System.out.println("2D RDF: Select the data within limits...");
		ArrayList<Integer> Z = new ArrayList<Integer>();
		for (int i = 0; i < Xa.length; i++) {
			// Removes 0,0 points which would otherwise skew the signal (as all particles
			// take the reference (0,0) position)
			if (Xa[i] == 0 && Ya[i] == 0) {
				continue;
			}
			// Only select the data within limits
			if (Xa[i] <= Rmax && Xa[i] >= -Rmax && Ya[i] <= Rmax && Ya[i] >= -Rmax) {
				Z.add(i);
			}
		}
		// Xa2 and Ya2 correspond to Xa and Ya in the matlab code
		float[] Xa2 = new float[Z.size()];
		float[] Ya2 = new float[Z.size()];
		for (int i = 0; i < Z.size(); i++) {
			Xa2[i] = Xa[Z.get(i)];
			Ya2[i] = Ya[Z.get(i)];
		}


		// Converts X,Y data from cartesian to polar- angle in radians
		double[] T_rad = new double[Xa2.length];
		double[] T = new double[Xa2.length];
		double[] R = new double[Xa2.length];

		for (int i = 0; i < Xa2.length; i++) {
			T_rad[i] = Math.atan2(Ya2[i], Xa2[i]);
			T[i] = T_rad[i];
			R[i] = Math.sqrt(Xa2[i] * Xa2[i] + Ya2[i] * Ya2[i]);
		}

		// Theta spacing was set in degrees, converted here to radians
		double Theta_Spacing_rad = spacing * Math.PI / 180;
		double[] Ts = createSpacedArray(0, Math.PI, Theta_Spacing_rad);

		// Pre-sets the variables for speed
		double[] Theta = new double[Ts.length - 1];
		double[] Tm = new double[Ts.length - 1];

		// Tm2 is set to plot 360 degrees based on rotational symmetry
		double[] Tm2 = new double[Ts.length - 1];
		System.arraycopy(Tm, 0, Tm2, 0, Tm.length);

		// Limits the data to within the Radial limits specified by user
		ArrayList<Double> T_temp = new ArrayList<>();
		for (int i = 0; i < R.length; i++) {
			if (R[i] >= Rmin && R[i] <= Rmax && R[i] > 0) {
				T_temp.add(T[i]);
			}
		}

		// Determines the number of points within each bin
		for (int i = 0; i < Ts.length - 1; i++) {
			for (int j = 0; j < T_temp.size(); j++) {
				if (T_temp.get(j) >= Ts[i] && T_temp.get(j) < Ts[i + 1]) {
					Theta[i] += 1;
				}
			}
			Tm[i] = (Ts[i + 1] + Ts[i]) / 2;
			Tm2[i] = Tm[i] + Math.PI;
		}

		// Connects the ends of the data for clarity in the polar plot
		double[] Tm_final = new double[Tm.length * 2];
		System.arraycopy(Tm, 0, Tm_final, 0, Tm.length);
		System.arraycopy(Tm2, 0, Tm_final, Tm.length, Tm2.length);

		double[] Theta_final = new double[Theta.length * 2];
		System.arraycopy(Theta, 0, Theta_final, 0, Theta.length);
		System.arraycopy(Theta, 0, Theta_final, Theta.length, Theta.length);

		// normalises the data
		double Theta_m = Arrays.stream(Theta_final).average().getAsDouble();
		double[] Theta_norm = new double[Theta_final.length];
		for (int i = 0; i < Theta_final.length; i ++) {
			Theta_norm[i]=Theta_final[i]/Theta_m;
		}

		// Plotting
		// First convert Tm to degree
		double[] Tm_degree = new double[Tm_final.length];
		for (int i =0; i < Tm_final.length; i++){
			Tm_degree[i] = Tm_final[i]/Math.PI*180;
		}
		XYSeries series = new XYSeries("Series 1");
		for (int i = 0; i < Theta_norm.length;i++) {
			series.add(Tm_degree[i],Theta_norm[i]);
		}
		XYSeriesCollection data = new XYSeriesCollection();
		data.addSeries(series);


		JFreeChart chart = ChartFactory.createPolarChart(
				null, data, false, false, false
		);
		PolarPlot plot = (PolarPlot) chart.getPlot();
		plot.setCounterClockwise(true);
		plot.setAxisLocation(PolarAxisLocation.EAST_BELOW);
		plot.setAngleOffset(0);
		plot.setBackgroundPaint(Color.white);
		Color gridColor = Color.lightGray;
		plot.setRadiusGridlinePaint(gridColor);
		plot.setAngleGridlinePaint(gridColor);
		DefaultPolarItemRenderer renderer = (DefaultPolarItemRenderer) plot.getRenderer();
		renderer.setShapesVisible(false);
		renderer.setSeriesFilled(0, false);
		renderer.setSeriesPaint(0,new Color((float) 0.0, (float) 0.4470, (float) 0.7410));
		ChartPanel chartPanel = new PolarChartPanel(chart);
		PlugInFrame frame = new PlugInFrame(null);
		frame.add(chartPanel);
		frame.pack();
		GUI.center(frame);
		frame.setVisible(true);
	}

	@SuppressWarnings("Duplicates")
	private float[] createSpacedArray(float min, float max, float spacing) {
		int count = (int) ((max - min)/spacing +1);
		float[] spacedArray = new float[count];

		for (int i = 0; i < count; i++) {
			spacedArray[i] = min + spacing * i;
		}
		return spacedArray;
	}
	@SuppressWarnings("Duplicates")
	private double[] createSpacedArray(double min, double max, double spacing) {
		int count = (int) ((max - min) / spacing + 1);
		double[] spacedArray = new double[count];

		for (int i = 0; i < count; i++) {
			spacedArray[i] = min + spacing * i;
		}
		return spacedArray;
	}
	@SuppressWarnings("Duplicates")
	private int[] floatArray2IntArray(float[] floatArray) {
		int[] intArray = new int[floatArray.length];
		for (int i = 0 ; i < floatArray.length; i++)
		{
			intArray[i] = (int) floatArray[i];
		}
		return intArray;
	}
}


