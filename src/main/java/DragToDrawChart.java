import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.ScheduledFuture;
import java.util.concurrent.TimeUnit;

public class DragToDrawChart extends JFrame {

    private final XYSeries sampleSeries;
    private final XYSeries functionSeries;
    private final List<Point2D.Double> drawnPoints = new ArrayList<>();
    private final XYPlot plot;
    private ScheduledFuture<?> future;
    private ScheduledExecutorService scheduledExecutorService;
    private double t1;
    private double t2;
    private double t3;
    private double t4;

    public DragToDrawChart() {
        super("Drag to Draw Points");

        sampleSeries = new XYSeries("Drawn Points");
        functionSeries = new XYSeries("polynomial_segmented_function");

        XYSeriesCollection dataset = new XYSeriesCollection(sampleSeries);
        dataset.addSeries(functionSeries);

        JFreeChart chart = ChartFactory.createXYLineChart(
                "Polyfit where segment-interval-times are variables",
                "X",
                "Y",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );


        plot = chart.getXYPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(true, false);
        plot.setRenderer(renderer);
        plot.getDomainAxis().setRange(0.0, 10.0);   // X axis
        plot.getRangeAxis().setRange(-10.0, 10.0);
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setMouseZoomable(false); // disable zoom

        chartPanel.addMouseMotionListener(new MouseAdapter() {
            @Override
            public void mouseDragged(MouseEvent e) {
                addPointFromMouse(e, chartPanel, plot);
            }
        });

        chartPanel.addMouseListener(new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                sampleSeries.clear();
                functionSeries.clear();

                drawnPoints.clear();
                addPointFromMouse(e, chartPanel, plot);
            }
        });


        setContentPane(chartPanel);
        setSize(800, 600);
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        setLocationRelativeTo(null);
    }

    private void addPointFromMouse(MouseEvent e, ChartPanel chartPanel, XYPlot plot) {
        Rectangle2D dataArea = chartPanel.getScreenDataArea();

        ValueAxis xAxis = plot.getDomainAxis();
        ValueAxis yAxis = plot.getRangeAxis();

        double x = xAxis.java2DToValue(
                e.getX(), dataArea, plot.getDomainAxisEdge());
        double y = yAxis.java2DToValue(
                e.getY(), dataArea, plot.getRangeAxisEdge());

        Point2D.Double point = new Point2D.Double(x, y);
        drawnPoints.add(point);
        sampleSeries.add(x, y);
        //functionSeries.add(x, y+1.0);

        scheduledExecutorService = Executors.newScheduledThreadPool(1);
        if (future != null) {
            future.cancel(false);
        }

        future =
                scheduledExecutorService.schedule(this::calculate, 1000L, TimeUnit.MILLISECONDS);

    }

    public List<Point2D.Double> getDrawnPoints() {
        return drawnPoints;
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            new DragToDrawChart().setVisible(true);
        });
    }

    private void calculate() {
        System.out.println("calculating");
        double[] p1_x = new double[]{1.02, 2.2, 3.123};
        double[] p1_y = new double[]{1.4, 5.765, 7.0};
        double[] p2_x = new double[]{11.0, 12.2, 12.0, 13.123, 13.2, 14.123};
        double[] p2_y = new double[]{11.4, 15.765, 11.4, 17.0, 15.765, 17.0};

        double[] p3_x = new double[]{18, 19, 20};
        double[] p3_y = new double[]{10, 15, 13};

        t1 = sampleSeries.getMinX();
        t2 = 3.123;
        t3 = 17;
        t4 = sampleSeries.getMaxX();

        double avgDelta= (sampleSeries.getMaxX()-sampleSeries.getMinX())/2;

        List<PolynomialFunction> result=null;
        double[] b_p1_x = new double[]{};
        double[] b_p1_y = new double[]{};
        double[] b_p2_x = new double[]{};
        double[] b_p2_y = new double[]{};
        double[] b_p3_x = new double[]{};
        double[] b_p3_y = new double[]{};
        double b_t2=0.0;
        double b_t3=0.0;
        List<PolynomialFunction> b_result=null;
        double best = 10000.0;

        for (int k=0;k<10000;k++){
            t2 =(Math.random()*avgDelta)+ t1;
            t3 =(Math.random()*avgDelta)+ t2;

            p1_x = sampleSeries.getItems().stream().filter(i->((XYDataItem)i).getX().doubleValue()>= t1 && ((XYDataItem)i).getX().doubleValue()< t2).mapToDouble(i->((XYDataItem)i).getX().doubleValue()).toArray();
            p2_x = sampleSeries.getItems().stream().filter(i->((XYDataItem)i).getX().doubleValue()>= t2 && ((XYDataItem)i).getX().doubleValue()< t3).mapToDouble(i->((XYDataItem)i).getX().doubleValue()).toArray();
            p3_x = sampleSeries.getItems().stream().filter(i->((XYDataItem)i).getX().doubleValue()>= t3 && ((XYDataItem)i).getX().doubleValue()< t4).mapToDouble(i->((XYDataItem)i).getX().doubleValue()).toArray();
            p1_y = sampleSeries.getItems().stream().filter(i->((XYDataItem)i).getX().doubleValue()>= t1 && ((XYDataItem)i).getX().doubleValue()< t2).mapToDouble(i->((XYDataItem)i).getY().doubleValue()).toArray();
            p2_y = sampleSeries.getItems().stream().filter(i->((XYDataItem)i).getX().doubleValue()>= t2 && ((XYDataItem)i).getX().doubleValue()< t3).mapToDouble(i->((XYDataItem)i).getY().doubleValue()).toArray();
            p3_y = sampleSeries.getItems().stream().filter(i->((XYDataItem)i).getX().doubleValue()>= t3 && ((XYDataItem)i).getX().doubleValue()< t4).mapToDouble(i->((XYDataItem)i).getY().doubleValue()).toArray();


            List<WeightedObservedPoint> samples1 = new ArrayList<>();
            for (int i = 0; i < p1_x.length; i++) {
                samples1.add(new WeightedObservedPoint(1.0, p1_x[i], p1_y[i]));
            }
            SegmentSampleData segment1 = new SegmentSampleData(t1, t2, samples1, 0, 0);

            List<WeightedObservedPoint> samples2 = new ArrayList<>();
            for (int i = 0; i < p2_x.length; i++) {
                samples2.add(new WeightedObservedPoint(1.0, p2_x[i], p2_y[i]));
            }
            SegmentSampleData segment2 = new SegmentSampleData(t2, t3, samples2, 30, 10);

            List<WeightedObservedPoint> samples3 = new ArrayList<>();
            for (int i = 0; i < p3_x.length; i++) {
                samples3.add(new WeightedObservedPoint(1.0, p3_x[i], p3_y[i]));
            }
            SegmentSampleData segment3 = new SegmentSampleData(t3, t4, samples3, 30, 10);

            List<SegmentSampleData> segments = List.of(segment1, segment2, segment3);

            PolyfitDto dto = new PolyfitDto(segments, 3);
            CustomPolyFitter customPolyFitter = new CustomPolyFitter(dto);
            //System.out.println(String.format("%f, %f, %f, %f", t1, t2, t3, t4));
            List<List<Double>> coeffs = customPolyFitter.calculateOptimalCoeffs();

            PolynomialFunction p1 = new PolynomialFunction(coeffs.get(0).stream().mapToDouble(Double::doubleValue).toArray());
            PolynomialFunction p2 = new PolynomialFunction(coeffs.get(1).stream().mapToDouble(Double::doubleValue).toArray());
            PolynomialFunction p3 = new PolynomialFunction(coeffs.get(2).stream().mapToDouble(Double::doubleValue).toArray());


            result = List.of(p1, p2, p3);
            if (CustomPolyFitter.best<best) {
                best=CustomPolyFitter.best;
                System.out.println(String.format("New best segmentation at: %f, %f, %f, %f , fit-value= %f", t1, t2, t3, t4,best));

                b_p1_x=p1_x;
                b_p1_y=p1_y;
                b_p2_x=p2_x;
                b_p2_y=p2_y;
                b_p3_x=p3_x;
                b_p3_y=p3_y;
                b_t2=t2;
                b_t3=t3;
                b_result=result;
            }
        }

        p1_x=b_p1_x;
        p1_y=b_p1_y;
        p2_x=b_p2_x;
        p2_y=b_p2_y;
        p3_x=b_p3_x;
        p3_y=b_p3_y;
        t2=b_t2;
        t3=b_t3;
        result=b_result;



        plotPolynomial(p1_x, result.get(0), plot, 4, "p1", Color.MAGENTA, t1, t2);
        plotPolynomial(p2_x, result.get(1), plot, 5, "p2", Color.ORANGE, t2, t3);
        plotPolynomial(p3_x, result.get(2), plot, 6, "p3", Color.RED, t3, t4);

        plotSamplePoints(plot, p1_x, p1_y, 1, Color.BLUE);
        plotSamplePoints(plot, p2_x, p2_y, 2, Color.YELLOW);
        plotSamplePoints(plot, p3_x, p3_y, 3, Color.PINK);


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