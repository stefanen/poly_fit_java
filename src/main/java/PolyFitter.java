import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
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
import java.util.List;

public class PolyFitter {

    public static final int DEGREE_TO_FIT = 3;

    public static void main(String[] args) {
        //simpleDemo1(1);
        //simpleDemo1(2);
        //simpleDemo1(3);
        //simpleDemo1(4);
        simpleDemo1(6);
    }


    private static void simpleDemo1(int demoCase) {
        double[] p1_x = {1.02, 2.2, 3.123};
        double[] p1_y = {1.4, 5.765, 7.0};
        double[] p2_x = {3.123, 11.0, 12.2, 12.0, 13.123, 13.2, 14.123};
        double[] p2_y = {7.0, 11.4, 15.765, 11.4, 17.0, 15.765, 17.0};


        XYPlot plot = new XYPlot();
        if (demoCase<5) {
            plotSamplePoints(plot, p1_x, p1_y, 1, Color.DARK_GRAY);
            plotSamplePoints(plot, p2_x, p2_y, 2, Color.BLUE);
        }
        /*
        PolynomialFunction p1Polynomial = PolyFitterUtil.fitSamplesDefault(p1_x, p1_y, DEGREE_TO_FIT);
        PolynomialFunction p2Polynomial = PolyFitterUtil.fitSamplesDefault(p2_x, p2_y, DEGREE_TO_FIT);
        */


        String title = "Polyfit";
        if (demoCase == 1) {
            PolynomialFunction p1Polynomial = PolyFitterUtil.fitSamplesDefault(p1_x, p1_y, DEGREE_TO_FIT);
            PolynomialFunction p2Polynomial = PolyFitterUtil.fitSamplesDefault(p2_x, p2_y, DEGREE_TO_FIT);
            plotPolynomial(p1_x, p1Polynomial, plot, 3, "p1", Color.LIGHT_GRAY);
            plotPolynomial(p2_x, p2Polynomial, plot, 4, "p2", Color.ORANGE);
            title += " reference";
        } else if (demoCase == 2) {
            Pair<PolynomialFunction, PolynomialFunction> fittedPolynomials = PolyFitterUtil.fitTwoAdjoiningPolynomialsSmoothly(p1_x, p1_y, p2_x, p2_y, DEGREE_TO_FIT, 10, 10);
            plotPolynomial(p1_x, fittedPolynomials.getFirst(), plot, 3, "p1", Color.LIGHT_GRAY);
            plotPolynomial(p2_x, fittedPolynomials.getSecond(), plot, 4, "p2", Color.ORANGE);
            title += " continuous and deriv-continuous";
        } else if (demoCase == 3) {
            Pair<PolynomialFunction, PolynomialFunction> fittedPolynomials = PolyFitterUtil.fitTwoAdjoiningPolynomialsSmoothly(p1_x, p1_y, p2_x, p2_y, DEGREE_TO_FIT, 1000, 0);
            plotPolynomial(p1_x, fittedPolynomials.getFirst(), plot, 3, "p1", Color.LIGHT_GRAY);
            plotPolynomial(p2_x, fittedPolynomials.getSecond(), plot, 4, "p2", Color.ORANGE);
            title += " continuous and not deriv-continuous";
        } else if (demoCase == 4) {
            Pair<PolynomialFunction, PolynomialFunction> fittedPolynomials = PolyFitterUtil.fitTwoAdjoiningPolynomialsSmoothly(p1_x, p1_y, p2_x, p2_y, DEGREE_TO_FIT, 0, 100);
            plotPolynomial(p1_x, fittedPolynomials.getFirst(), plot, 3, "p1", Color.LIGHT_GRAY);
            plotPolynomial(p2_x, fittedPolynomials.getSecond(), plot, 4, "p2", Color.ORANGE);
            title += " not-continuous but deriv-continuous";
        } else if (demoCase == 5) {
            p1_x = new double[]{1.02, 2.2, 2};
            p1_y = new double[]{1.4, 5.765, 4};
            p2_x = new double[]{11.0, 12.2, 12.0, 13.123, 13.2, 14.123};
            p2_y = new double[]{11.4, 15.765, 11.4, 17.0, 15.765, 17.0};

            double t1=1;
            double t2=3.123;
            double t3=17;

            Pair<PolynomialFunction, PolynomialFunction> fittedPolynomials = PolyFitterUtil.fitTwo(p1_x, p1_y, p2_x, p2_y, DEGREE_TO_FIT, 3, 3,t1,t2,t3);
            plotPolynomial(p1_x, fittedPolynomials.getFirst(), plot, 3, "p1", Color.LIGHT_GRAY,t1,t2);
            plotPolynomial(p2_x, fittedPolynomials.getSecond(), plot, 4, "p2", Color.ORANGE, t2,t3);
            title += " 2 segments joining";
            plotSamplePoints(plot, p1_x, p1_y, 1, Color.DARK_GRAY);
            plotSamplePoints(plot, p2_x, p2_y, 2, Color.BLUE);

        } else if (demoCase == 6) {
            p1_x = new double[]{1.02, 2.2, 3.123};
            p1_y = new double[]{1.4, 5.765, 7.0};
            p2_x = new double[]{11.0, 12.2, 12.0, 13.123, 13.2, 14.123};
            p2_y = new double[]{11.4, 15.765, 11.4, 17.0, 15.765, 17.0};

            double[] p3_x = new double[]{18, 19, 20};
            double[] p3_y = new double[]{10, 15, 13};

            double t1=1;
            double t2=3.123;
            double t3=17;
            double t4=21;



/*
            PolyFitterUtil.SegmentSampleData s1 = new PolyFitterUtil.SegmentSampleData(p1xSamples, p1ySamples, s1start, s1end);
            PolyFitterUtil.SegmentSampleData s2 = new PolyFitterUtil.SegmentSampleData(p2xSamples, p2ySamples, s1end, s2end);
            PolyFitterUtil.SegmentSampleData s3 = new PolyFitterUtil.SegmentSampleData(p3xSamples, p3ySamples, s2end, s3end);

            //SegmentSampleData s2 = new SegmentSampleData(p2xSamples,p2ySamples,p1xSamples[p1xSamples.length-1]), p2xSamples[p2xSamples.length-1]));
            List<List<Double>> coeffs2 = polyfitMultiple(List.of(s1), degree + 1, cWeight, dWeight);
            System.out.println(coeffs2);

            List<List<Double>> coeffs = polyfitMultiple(List.of(s1, s2, s3), degree + 1, cWeight, dWeight);
            PolynomialFunction p1 = new PolynomialFunction(coeffs.get(0).stream().mapToDouble(Double::doubleValue).toArray());
            PolynomialFunction p2 = new PolynomialFunction(coeffs.get(1).stream().mapToDouble(Double::doubleValue).toArray());
            PolynomialFunction p3 = new PolynomialFunction(coeffs.get(2).stream().mapToDouble(Double::doubleValue).toArray());
  */


            List<WeightedObservedPoint> samples1 = new ArrayList<>();
            for (int i=0; i<p1_x.length;i++) {
                samples1.add(new WeightedObservedPoint(1.0,p1_x[i],p1_y[i]));
            }
            SegmentSampleData segment1 = new SegmentSampleData(t1, t2, samples1, 0.3,3);

            List<WeightedObservedPoint> samples2 = new ArrayList<>();
            for (int i=0; i<p2_x.length;i++) {
                samples2.add(new WeightedObservedPoint(1.0,p2_x[i],p2_y[i]));
            }
            SegmentSampleData segment2 = new SegmentSampleData(t2, t3, samples2, 0.3,3);

            List<WeightedObservedPoint> samples3 = new ArrayList<>();
            for (int i=0; i<p3_x.length;i++) {
                samples3.add(new WeightedObservedPoint(1.0,p3_x[i],p3_y[i]));
            }
            SegmentSampleData segment3 = new SegmentSampleData(t3, t4, samples3, 3,3);

            List<SegmentSampleData> segments = List.of(segment1,segment2,segment3);

            PolyfitDto dto = new PolyfitDto(segments,DEGREE_TO_FIT);
            CustomPolyFitter customPolyFitter = new CustomPolyFitter(dto);

            List<List<Double>> coeffs = customPolyFitter.calculateOptimalCoeffs();
            PolynomialFunction p1 = new PolynomialFunction(coeffs.get(0).stream().mapToDouble(Double::doubleValue).toArray());
            PolynomialFunction p2 = new PolynomialFunction(coeffs.get(1).stream().mapToDouble(Double::doubleValue).toArray());
            PolynomialFunction p3 = new PolynomialFunction(coeffs.get(2).stream().mapToDouble(Double::doubleValue).toArray());

            List<PolynomialFunction> result = List.of(p1,p2,p3);





            plotPolynomial(p1_x, result.get(0), plot, 4, "p1", Color.LIGHT_GRAY, t1,t2);
            plotPolynomial(p2_x, result.get(1), plot, 5, "p2", Color.ORANGE, t2,t3);
            plotPolynomial(p3_x, result.get(2), plot, 6, "p3", Color.RED, t3,t4);

            plotSamplePoints(plot, p1_x, p1_y, 1, Color.DARK_GRAY);
            plotSamplePoints(plot, p2_x, p2_y, 2, Color.BLUE);
            plotSamplePoints(plot, p3_x, p3_y, 3, Color.DARK_GRAY);

            title += " 3 segments joining";
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



    private static void plotPolynomial(double[] x, PolynomialFunction polynomial, XYPlot plot, int id, String name, Color color, double start, double end) {
        java.util.List<Double> xValues = new ArrayList<>();
        java.util.List<Double> yValues = new ArrayList<>();

        XYSeries series2 = new XYSeries(name);
        for (double i = start; i < end; i += 0.01) {

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
