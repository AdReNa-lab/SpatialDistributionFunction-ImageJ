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
//
//VIB Protocol
//Copyright (C) 2017  J. Schindelin, Benjamin Schmid, M. Longair
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
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.MultiLineLabel;
import ij.measure.ResultsTable;
import ij.plugin.LutLoader;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import utils.GUIHelper;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;

@SuppressWarnings("Duplicates")
public class SpatialDistributionFunction_2D implements PlugInFilter {

    ImagePlus imp;
    private double XMaxDefault;
    private double YMaxDefault;
    private double SpacingDefault;

    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        return DOES_ALL + NO_UNDO;
    }

    public void run(ImageProcessor ip) {
        ResultsTable coordinateTable = Analyzer.getResultsTable();
        // Terminate here if table is not available. otherwise error will occur
        int xColumnIndex = coordinateTable.getColumnIndex("X");
        int yColumnIndex = coordinateTable.getColumnIndex("Y");
        // Read in positional information
        float[] xValues = coordinateTable.getColumn(xColumnIndex);
        float[] yValues = coordinateTable.getColumn(yColumnIndex);
        // Read from dialog
        double Xmax;
        double Ymax;
        float spacing;
        GenericDialog gd = new GenericDialog("2D RDF", IJ.getInstance());
        gd.addNumericField("X max (pixels):", XMaxDefault, 2);
        gd.addNumericField("Y max (pixels):", YMaxDefault, 2);
        gd.addNumericField("Spacing", SpacingDefault, 2);
        gd.addMessage("This plugin analyses colloidal ordering via the SDF and plots the resulting analysis in Cartesian and polar coordinates. \n(see scientific publication in Langmuir:  https://doi.org/10.1021/acs.langmuir.9b02877)");
        MultiLineLabel text = (MultiLineLabel) gd.getMessage();
        GUIHelper.addHyperLinkListener(text, "https://doi.org/10.1021/acs.langmuir.9b02877");
        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        }
        XMaxDefault = gd.getNextNumber();
        YMaxDefault = gd.getNextNumber();
        SpacingDefault = gd.getNextNumber();
        Xmax = XMaxDefault;
        Ymax = YMaxDefault;
        spacing = (float) SpacingDefault;
        if (gd.invalidNumber()) {
            IJ.showMessage("Error", "Invalid input Number");
            return;
        }

        float[] Xd = xValues.clone();
        float[] Yd = yValues.clone();
        // Sets up the edges for the 2D bins
        float[] Xs1 = createSpacedArray(0, (float) Xmax, spacing);
        float[] Xs2 = createSpacedArray(-spacing, (float) -Xmax, -spacing);
        float[] Xs = new float[Xs1.length + Xs2.length];
        for (int i = 0; i < Xs2.length; i++) {
            Xs[i] = Xs2[Xs2.length - 1 - i];
        }
        for (int i = 0; i < Xs1.length; i++) {
            Xs[i + Xs2.length] = Xs1[i];
        }

        float[] Ys1 = createSpacedArray(0, (float) Ymax, spacing);
        float[] Ys2 = createSpacedArray(-spacing, (float) -Ymax, -spacing);
        float[] Ys = new float[Ys1.length + Ys2.length];
        for (int i = 0; i < Ys2.length; i++) {
            Ys[i] = Ys2[Ys2.length - 1 - i];
        }
        for (int i = 0; i < Ys1.length; i++) {
            Ys[i + Ys2.length] = Ys1[i];
        }

        // Pre-sets size of matrix of particle density (2D-histogram)
        float[] H = new float[(Xs.length - 1) * (Ys.length - 1)];

        //Determines the dimensions of the sample
        float max_Xd = Xd[0];
        float min_Xd = Xd[0];
        for (int i = 1; i < Xd.length; i++) {
            if (Xd[i] > max_Xd) {
                max_Xd = Xd[i];
            }
            if (Xd[i] < min_Xd) {
                min_Xd = Xd[i];
            }
        }
        float max_Yd = Yd[0];
        float min_Yd = Yd[0];
        for (int i = 1; i < Yd.length; i++) {
            if (Yd[i] > max_Yd) {
                max_Yd = Yd[i];
            }
            if (Yd[i] < min_Yd) {
                min_Yd = Yd[i];
            }
        }
        float x_dim = max_Xd - min_Xd;
        float y_dim = max_Yd - min_Yd;
        // Sets the lower left position of the rectangular sample area
        float ref_pos1 = min_Xd;
        float ref_pos2 = min_Yd;
        // Presets size of matrix for normalisation
        float[] norm_mat = new float[(Ys.length - 1) * (Xs.length - 1)];
        // Calculate the average density overall
        float avg_density_overall = Xs.length / x_dim / y_dim;
        // Calculates the average density expected in a given bin
        float avg_density = avg_density_overall * spacing * spacing;

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
            for (int k = 0; k < Yd.length; k++) {
                Xa[i * Yd.length + k] = X_temp[k];
                Ya[i * Yd.length + k] = Y_temp[k];
            }

            // Shifts the rectangle (the sample area)
            float ref_pos_temp1 = ref_pos1 - Xc;
            float ref_pos_temp2 = ref_pos2 - Yc;

            //determines the limits of the shifted rectangle
            float rect_x_min = ref_pos_temp1;
            float rect_x_max = ref_pos_temp1 + x_dim;
            float rect_y_min = ref_pos_temp2;
            float rect_y_max = ref_pos_temp2 + y_dim;

            // Determines the area of the rectangle (comprising the sample area) that
            // overlaps with each cell in the grid produced for the SDF
            for (int j = 0; j < Xs.length - 1; j++) {
                float rx;
                float ry;
                if (rect_x_min >= Xs[j] && rect_x_min <= Xs[j + 1]) {
                    rx = Xs[j + 1] - rect_x_min;
                } else if (rect_x_min <= Xs[j] && rect_x_max >= Xs[j + 1]) {
                    rx = Xs[j + 1] - Xs[j];
                } else if (rect_x_max >= Xs[j] && rect_x_max <= Xs[j + 1]) {
                    rx = rect_x_max - Xs[j];
                } else {
                    rx = 0;
                }
                for (int k = 0; k < Ys.length - 1; k++) {
                    if (rect_y_min >= Ys[k] && rect_y_min <= Ys[k + 1]) {
                        ry = Ys[k + 1] - rect_y_min;
                    } else if (rect_y_min <= Ys[k] && rect_y_max >= Ys[k + 1]) {
                        ry = Ys[k + 1] - Ys[k];
                    } else if (rect_y_max >= Ys[k] && rect_y_max <= Ys[k + 1]) {
                        ry = rect_y_max - Ys[k];
                    } else {
                        ry = 0;
                    }
                    float ra = rx * ry / spacing / spacing;
                    norm_mat[k * (Xs.length - 1) + j] += ra;
                }

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
            if (Xa[i] <= Xmax && Xa[i] >= -Xmax && Ya[i] <= Ymax && Ya[i] >= -Ymax) {
                Z.add(i);
            }
        }
        float[] Xa2 = new float[Z.size()];
        float[] Ya2 = new float[Z.size()];
        for (int i = 0; i < Z.size(); i++) {
            Xa2[i] = Xa[Z.get(i)];
            Ya2[i] = Ya[Z.get(i)];
        }

//		// Sets up the edges for the 2D bins
//		System.out.println("2D RDF: Sets up the edges for the 2D bins...");
//		float[] Xs = createSpacedArray((float) -Xmax, (float) Xmax,spacing);
//		float[] Ys = createSpacedArray((float) -Ymax, (float) Ymax,spacing);


        for (int i = 0; i < Xs.length - 1; i++) {
            ArrayList<Integer> X3 = new ArrayList<Integer>();
            for (int q = 0; q < Xa2.length; q++) {
                if (Xa2[q] >= Xs[i] && Xa2[q] < Xs[i + 1]) {
                    X3.add(q);
                }
            }

            for (int j = 0; j < Ys.length - 1; j++) {
                ArrayList<Integer> Z3 = new ArrayList<Integer>();
                for (int w = 0; w < Ya2.length; w++) {
                    if (Ya2[w] < Ys[j + 1] && Ya2[w] >= Ys[j] && X3.contains(w)) {
                        Z3.add(w);
                    }
                }
                H[i + j * (Xs.length - 1)] = Z3.size() / spacing / spacing;
            }
        }

        // Calculates the mid-point of each bin for plotting purposes
        System.out.println("2D RDF: Calculates the mid-point of each bin...");
        float[] Xm = new float[Xs.length - 1];
        for (int i = 0; i < Xs.length - 1; i++) {
            Xm[i] = (Xs[i] + Xs[i + 1]) / 2;
        }
        float[] Ym = new float[Ys.length - 1];
        for (int i = 0; i < Ys.length - 1; i++) {
            Ym[i] = (Ys[i] + Ys[i + 1]) / 2;
        }

        // normalizes data
        System.out.println("2D RDF: Normalizes data...");
        float[] Hm = new float[H.length];
//		float m;
//		float sumOfH = 0;
//		for (int i = 0; i < H.length; i ++) {
//			sumOfH += H[i];
//		}
//		m = sumOfH/H.length;
//		for (int i = 0; i < H.length; i ++ ){
//			Hm[i]=H[i]/m;
//		}
        for (int i = 0; i < norm_mat.length; i++) {
            if (norm_mat[i] == 0) {
                norm_mat[i] = 1;
                Hm[i] = 0;
            } else {
                Hm[i] = H[i] / (norm_mat[i] * avg_density);
            }
        }

        System.out.println("2D RDF: Calculation finished");

        // Convert to byte image
        float HmMin = Hm[0];
        float HmMax = Hm[0];
        for (int i = 0; i < Hm.length; i++) {
            if (HmMin > Hm[i]) {
                HmMin = Hm[i];
            }
            if (HmMax < Hm[i]) {
                HmMax = Hm[i];
            }
        }
        byte[] HmByte = new byte[Hm.length];
        for (int i = 0; i < Hm.length; i++) {
            HmByte[i] = (byte) ((Hm[i] - HmMin) / (HmMax - HmMin) * 255);
        }


//		float[] testvals = new float[900];
//
//		for (int x = 0; x< 900; x++){
//			testvals[x] = (float) x / 300;
//		}
        ByteProcessor fp = new ByteProcessor(Xs.length - 1, Ys.length - 1);
        fp.setPixels(HmByte);


        ImagePlus testPlot;

        if (Xmax * 2 / spacing < 150 && Ymax * 2 / spacing < 150) {
            testPlot = new ImagePlus("SDF", fp.resize(300));
        } else {
            testPlot = new ImagePlus("SDF", fp);
        }

//		testPlot.show();

        ImagePlus newPlot = addAxes(testPlot, "X", (double) -Xmax, (double) Xmax,
                "Y", (double) -Ymax, (double) Ymax, HmMin, HmMax);

        newPlot.show();
    }

//	public ImagePlus addFrame(
//			String title,
//			ImagePlus baseRdf,
//			String xLabel, double xmin, double xmax,
//			String yLabel, double ymin, double ymax) {
//
//	}

    public ImagePlus addAxes(
            ImagePlus baseRdf,
            String xLabel, double xmin, double xmax,
            String yLabel, double ymin, double ymax,
            double barMin, double barMax) {

        int tickSize = 5;
        int tickMargin = 10;

        int leftBorder = 100;
        int rightBorder = 180;
        int topBorder = 60;
        int bottomBorder = 100;


        int oldWidth = baseRdf.getWidth();
        int oldHeight = baseRdf.getHeight();
        ByteProcessor oldFP = (ByteProcessor) baseRdf.getProcessor();
        double oldMin = (double) oldFP.getMin();
        double oldMax = (double) oldFP.getMax();
        byte[] oldFloats = (byte[]) oldFP.getPixels();

        int newWidth = oldWidth + leftBorder + rightBorder;
        int newHeight = oldHeight + topBorder + bottomBorder;
        byte[] newFloats = new byte[newWidth * newHeight];
        for (int i = 0; i < newFloats.length; ++i)
            newFloats[i] = (byte) oldMax;

        for (int y = 0; y < oldHeight; ++y) {
            for (int x = 0; x < oldWidth; ++x) {
                newFloats[(y + topBorder) * newWidth + (x + leftBorder)] =
                        oldFloats[y * oldWidth + x];
            }
        }

        ByteProcessor newFP = new ByteProcessor(newWidth, newHeight);
        newFP.setPixels(newFloats);
        newFP.setMinAndMax(oldMin, oldMax);

        newFP.setValue(oldMin);

        // Draw ticks:
        // Ticks on x axis
        newFP.drawLine(
                leftBorder,
                topBorder + oldHeight,
                leftBorder,
                topBorder + oldHeight + tickSize);
        newFP.drawLine(
                leftBorder + oldWidth - 1,
                topBorder + oldHeight,
                leftBorder + oldWidth - 1,
                topBorder + oldHeight + tickSize);
        newFP.drawLine(
                (int) (leftBorder + oldWidth / 6 - 1),
                topBorder + oldHeight,
                (int) (leftBorder + oldWidth / 6 - 1),
                topBorder + oldHeight + tickSize);
        newFP.drawLine(
                (int) (leftBorder + oldWidth / 3 - 1),
                topBorder + oldHeight,
                (int) (leftBorder + oldWidth / 3 - 1),
                topBorder + oldHeight + tickSize);
        newFP.drawLine(
                (int) (leftBorder + oldWidth / 2 - 1),
                topBorder + oldHeight,
                (int) (leftBorder + oldWidth / 2 - 1),
                topBorder + oldHeight + tickSize);
        newFP.drawLine(
                (int) (leftBorder + oldWidth * 2 / 3 - 1),
                topBorder + oldHeight,
                (int) (leftBorder + oldWidth * 2 / 3 - 1),
                topBorder + oldHeight + tickSize);
        newFP.drawLine(
                (int) (leftBorder + oldWidth * 5 / 6 - 1),
                topBorder + oldHeight,
                (int) (leftBorder + oldWidth * 5 / 6 - 1),
                topBorder + oldHeight + tickSize);


        // Ticks on y axis
        newFP.drawLine(
                leftBorder - 1,
                topBorder,
                (leftBorder - 1) - tickSize,
                topBorder);
        newFP.drawLine(
                leftBorder - 1,
                topBorder + oldHeight - 1,
                (leftBorder - 1) - tickSize,
                topBorder + oldHeight - 1);
        newFP.drawLine(
                leftBorder - 1,
                (int) (topBorder + oldHeight / 6 - 1),
                (leftBorder - 1) - tickSize,
                (int) (topBorder + oldHeight / 6 - 1));
        newFP.drawLine(
                leftBorder - 1,
                (int) (topBorder + oldHeight / 3 - 1),
                (leftBorder - 1) - tickSize,
                (int) (topBorder + oldHeight / 3 - 1));
        newFP.drawLine(
                leftBorder - 1,
                (int) (topBorder + oldHeight / 2 - 1),
                (leftBorder - 1) - tickSize,
                (int) (topBorder + oldHeight / 2 - 1));
        newFP.drawLine(
                leftBorder - 1,
                (int) (topBorder + oldHeight * 2 / 3 - 1),
                (leftBorder - 1) - tickSize,
                (int) (topBorder + oldHeight * 2 / 3 - 1));
        newFP.drawLine(
                leftBorder - 1,
                (int) (topBorder + oldHeight * 5 / 6 - 1),
                (leftBorder - 1) - tickSize,
                (int) (topBorder + oldHeight * 5 / 6 - 1));

        ImagePlus newImagePlus = new ImagePlus(
                "SDF 2D",
                newFP);


        newImagePlus.getProcessor().setColorModel(baseRdf.getProcessor().getColorModel());

//		serifFont = true;
//		String fontName = serifFont ? "Serif" : "SanSerif";
//		int fontType = false ? Font.BOLD : Font.PLAIN;

        // This is font for the tick numbers
        Font font = new Font("Arial", Font.PLAIN, 14);

        newImagePlus.show();
        ImageCanvas ic = newImagePlus.getCanvas();
        FontMetrics fm = ic.getFontMetrics(font);

        newFP.setFont(font);
        newFP.setAntialiasedText(false);

        String sXmin = "" + xmin;
        String sXmax = "" + xmax;
        String sYmin = "" + ymin;
        String sYmax = "" + ymax;

        newFP.drawString(
                sXmin,
                leftBorder - (fm.stringWidth(sXmin) / 2),
                topBorder + oldHeight + tickSize + tickMargin + fm.getHeight());
        newFP.drawString(
                sXmax,
                leftBorder + oldWidth - (fm.stringWidth(sXmax) / 2),
                topBorder + oldHeight + tickSize + tickMargin + fm.getHeight());
        newFP.drawString(
                "" + (xmin + xmax) / 2,
                leftBorder + oldWidth / 2 - (fm.stringWidth(sXmax) / 2),
                topBorder + oldHeight + tickSize + tickMargin + fm.getHeight());
        newFP.drawString(
                sYmin,
                leftBorder - tickMargin - fm.stringWidth(sYmin) - tickSize,
                topBorder + oldHeight + fm.getHeight() / 2);
        newFP.drawString(
                sYmax,
                leftBorder - tickMargin - fm.stringWidth(sYmax) - tickSize,
                topBorder + fm.getHeight() / 2);
        newFP.drawString(
                "" + (ymax + ymin) / 2,
                leftBorder - tickMargin - fm.stringWidth(sYmin) - tickSize,
                topBorder + oldHeight / 2 + fm.getHeight() / 2);


        Font axisLabelFont = new Font("Arial", Font.PLAIN, 18);

        newImagePlus.show();
        FontMetrics axisLabelFm = ic.getFontMetrics(axisLabelFont);

        newFP.setFont(axisLabelFont);
        newFP.setAntialiasedText(false);

        newFP.drawString(
                xLabel,
                leftBorder + oldWidth / 2 - axisLabelFm.stringWidth(xLabel) / 2,
                topBorder + oldHeight + tickSize + 2 * tickMargin + 2 * axisLabelFm.getHeight());

        /* Draw a similar label in a new FloatProcessor and copy
         * it over. */

        int labelWidth = axisLabelFm.stringWidth(yLabel);
        int labelHeight = axisLabelFm.getHeight();

        ByteProcessor fpToRotate = new ByteProcessor(labelWidth, labelHeight);
        byte[] labelFloats = new byte[labelWidth * labelHeight];
        for (int i = 0; i < labelFloats.length; ++i)
            labelFloats[i] = (byte) oldMax;
        fpToRotate.setFont(axisLabelFont);
        fpToRotate.setPixels(labelFloats);
        fpToRotate.setValue(oldMin);
        fpToRotate.setMinAndMax(oldMin, oldMax);
        fpToRotate.drawString(yLabel, 0, labelHeight);

        int yLabelTopLeftX = leftBorder - tickSize - 2 * tickMargin - labelHeight * 2 - fm.getHeight();
        int yLabelTopLeftY = topBorder + (oldHeight / 2) - (labelWidth / 2);

        for (int y = 0; y < labelHeight; ++y)
            for (int x = 0; x < labelWidth; ++x) {
                int newX = yLabelTopLeftX + y;
                int newY = yLabelTopLeftY + labelWidth - x;
                newFloats[newY * newWidth + newX] = labelFloats[y * labelWidth + x];
            }

        /* Now draw a bar at the side showing the value range. */

        int barWidth = 30;
        int barHeight = (oldHeight * 2) / 3;

        int barTopLeftX = leftBorder + oldWidth + 40;
        int barTopLeftY = topBorder + (oldHeight - barHeight) / 2;

        newFP.drawRect(barTopLeftX, barTopLeftY, barWidth + 2, barHeight + 2);

        for (int barOffset = 0; barOffset < barHeight; ++barOffset) {
            int barLineX1 = barTopLeftX + 1;
            int barLineX2 = barTopLeftX + barWidth;
            int barLineY = barTopLeftY + 1 + (barHeight - (barOffset + 1));
            double value = ((barOffset * (oldMax - oldMin)) / (barHeight - 1) + oldMin);
//			byte value= (byte) ((barOffset*255)/(barHeight-1));
            newFP.setValue(value);
            newFP.drawLine(barLineX1, barLineY, barLineX2, barLineY);
        }

        /* Now add some tick marks to the bar */
        // use the font size of tick numbers
        newFP.setFont(font);

        DecimalFormat df = new DecimalFormat("0.0");
        newFP.setValue(oldMin);
        newFP.drawLine(
                barTopLeftX + barWidth + 2,
                barTopLeftY,
                barTopLeftX + barWidth + 2 + tickSize,
                barTopLeftY);
        newFP.drawString(
                "" + df.format(barMax),
                barTopLeftX + barWidth + 2 + tickSize + tickMargin,
                barTopLeftY + fm.getHeight() / 2
        );
        newFP.drawLine(
                barTopLeftX + barWidth + 2,
                barTopLeftY + barHeight + 1,
                barTopLeftX + barWidth + 2 + tickSize,
                barTopLeftY + barHeight + 1);
        newFP.drawString(
                "" + df.format(barMin),
                barTopLeftX + barWidth + 2 + tickSize + tickMargin,
                barTopLeftY + barHeight + fm.getHeight() / 2
        );

        /* Now just draw the title */

//		fontType = Font.BOLD;
//		Font titleFont=new Font(fontName, fontType, titleSize);

//		FontMetrics titleFM=ic.getFontMetrics(font);

//		newFP.setFont(titleFont);
//		newFP.drawString(
//				title,
//				newWidth / 2 - titleFM.stringWidth(title) / 2,
//				topBorder / 2 + titleFM.getHeight() / 2 );

        newImagePlus.updateAndRepaintWindow();


        newImagePlus.setLut(LutLoader.openLut("luts/Red Hot.lut"));

        return newImagePlus;
    }


    private float[] createSpacedArray(float min, float max, float spacing) {
        int count = (int) ((max - min) / spacing + 1);
        float[] spacedArray = new float[count];

        for (int i = 0; i < count; i++) {
            spacedArray[i] = min + spacing * i;
        }
        return spacedArray;
    }


    private int[] floatArray2IntArray(float[] floatArray) {
        int[] intArray = new int[floatArray.length];
        for (int i = 0; i < floatArray.length; i++) {
            intArray[i] = (int) floatArray[i];
        }
        return intArray;
    }


}


