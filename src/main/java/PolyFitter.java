import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.util.Pair;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;

public class PolyFitter {

    public static final int DEGREE_TO_FIT = 3;

    public static void main(String[] args) {
        //simpleDemo1(1);
        simpleDemo1(2);
        simpleDemo1(3);
        simpleDemo1(4);
    }


    private static void simpleDemo1(int demoCase) {
        double[] p1_x = {1.02, 2.2, 3.123};
        double[] p1_y = {1.4, 5.765, 7.0};
        double[] p2_x = {3.123, 11.0, 12.2, 12.0, 13.123, 13.2, 14.123};
        double[] p2_y = {7.0, 11.4, 15.765, 11.4, 17.0, 15.765, 17.0};


        XYPlot plot = new XYPlot();
        plotSamplePoints(plot, p1_x, p1_y, 1, Color.DARK_GRAY);
        plotSamplePoints(plot, p2_x, p2_y, 2, Color.BLUE);

        /*
        PolynomialFunction p1Polynomial = PolyFitterUtil.fitSamplesDefault(p1_x, p1_y, DEGREE_TO_FIT);
        PolynomialFunction p2Polynomial = PolyFitterUtil.fitSamplesDefault(p2_x, p2_y, DEGREE_TO_FIT);
        */

        String title="Polyfit";
        if (demoCase==1) {
            PolynomialFunction p1Polynomial = PolyFitterUtil.fitSamplesDefault(p1_x, p1_y, DEGREE_TO_FIT);
            PolynomialFunction p2Polynomial = PolyFitterUtil.fitSamplesDefault(p2_x, p2_y, DEGREE_TO_FIT);
            plotPolynomial(p1_x, p1Polynomial, plot, 3, "p1", Color.LIGHT_GRAY);
            plotPolynomial(p2_x, p2Polynomial, plot, 4, "p2", Color.ORANGE);
            title+=" reference";
        } else if (demoCase==2) {
            Pair<PolynomialFunction, PolynomialFunction> fittedPolynomials = PolyFitterUtil.fitTwoAdjoiningPolynomialsSmoothly(p1_x, p1_y, p2_x, p2_y, DEGREE_TO_FIT,10,10);
            plotPolynomial(p1_x, fittedPolynomials.getFirst(), plot, 3, "p1", Color.LIGHT_GRAY);
            plotPolynomial(p2_x, fittedPolynomials.getSecond(), plot, 4, "p2", Color.ORANGE);
            title+=" continuous and deriv-continuous";
        } else if (demoCase==3) {
            Pair<PolynomialFunction, PolynomialFunction> fittedPolynomials = PolyFitterUtil.fitTwoAdjoiningPolynomialsSmoothly(p1_x, p1_y, p2_x, p2_y, DEGREE_TO_FIT,1000,0);
            plotPolynomial(p1_x, fittedPolynomials.getFirst(), plot, 3, "p1", Color.LIGHT_GRAY);
            plotPolynomial(p2_x, fittedPolynomials.getSecond(), plot, 4, "p2", Color.ORANGE);
            title+=" continuous and not deriv-continuous";
        } else if (demoCase==4) {
            Pair<PolynomialFunction, PolynomialFunction> fittedPolynomials = PolyFitterUtil.fitTwoAdjoiningPolynomialsSmoothly(p1_x, p1_y, p2_x, p2_y, DEGREE_TO_FIT,0,100);
            plotPolynomial(p1_x, fittedPolynomials.getFirst(), plot, 3, "p1", Color.LIGHT_GRAY);
            plotPolynomial(p2_x, fittedPolynomials.getSecond(), plot, 4, "p2", Color.ORANGE);
            title+=" not-continuous but deriv-continuous";
        }






        JFreeChart chart = new JFreeChart(
                title,
                null,
                plot,
                true);


        JFrame frame = new JFrame("Polyfit experiment");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.pack();
        frame.setVisible(true);
    }

    private static void plotPolynomial(double[] x, PolynomialFunction polynomial, XYPlot plot, int id, String name, Color color) {
        java.util.List<Double> xValues = new ArrayList<>();
        java.util.List<Double> yValues = new ArrayList<>();

        XYSeries series2 = new XYSeries(name);
        for (double i = x[0]; i < x[x.length - 1]; i += 0.01) {

            series2.add(i, polynomial.value(i));
            xValues.add(i);
            yValues.add(polynomial.value(i));
        }

        XYSeriesCollection dataset2 = new XYSeriesCollection(series2);
        XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer(true, true);
        renderer2.setSeriesPaint(0, color);
        plot.setDataset(id, new XYSeriesCollection(series2));
        plot.setRenderer(id, renderer2);


        plot.setDatasetRenderingOrder(DatasetRenderingOrder.REVERSE);
    }

    private static XYPlot plotSamplePoints(XYPlot plot, double[] x, double[] y, int id, Color color) {
        XYSeries series = new XYSeries("Samples from p_" + id);
        for (int i = 0; i < x.length; i++) {
            series.add(x[i], y[i]);
        }
        XYSeriesCollection dataset = new XYSeriesCollection(series);


        NumberAxis xAxis = new NumberAxis("X-Axis");
        NumberAxis yAxis = new NumberAxis("Y-Axis");
        plot.setDomainAxis(xAxis);
        plot.setRangeAxis(yAxis);

        // Create renderers for different datasets
        XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer(false, true); // Only shapes (scatter)
        renderer1.setSeriesPaint(0, color);
        renderer1.setSeriesPaint(1, color);
        plot.setDataset(id, dataset);
        plot.setRenderer(id, renderer1);
        return plot;
    }

}
