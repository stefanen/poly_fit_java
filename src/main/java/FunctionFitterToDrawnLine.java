import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
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
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.ScheduledFuture;
import java.util.concurrent.TimeUnit;

public class FunctionFitterToDrawnLine extends JFrame {

    private final XYSeries sampleSeries;
    private final List<Point2D.Double> drawnPoints = new ArrayList<>();
    private final XYPlot plot;
    private ScheduledFuture<?> calculationFuture;

    public FunctionFitterToDrawnLine() {
        super("Drag to Draw Points");

        sampleSeries = new XYSeries("Drawn Points");

        XYSeriesCollection dataset = new XYSeriesCollection(sampleSeries);

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

        ScheduledExecutorService scheduledExecutorService = Executors.newScheduledThreadPool(1);
        if (calculationFuture != null) {
            calculationFuture.cancel(false);
        }

        calculationFuture =
                scheduledExecutorService.schedule(this::calculate, 1000L, TimeUnit.MILLISECONDS);

    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            new FunctionFitterToDrawnLine().setVisible(true);
        });
    }

    record PartitionedSample(double startTime, double endTime, List<Point2D.Double> points) {

    }


    private PartitionedSample getPartition(double startTime, double endTime, XYSeries allSamples) {
        List<XYDataItem> items = sampleSeries.getItems().stream().map(i->(XYDataItem)i).toList();

        return new PartitionedSample(startTime,endTime,items.stream()
                .filter((XYDataItem i)->i.getX().doubleValue()>= startTime && i.getX().doubleValue()< endTime)
                .map((XYDataItem i)->new Point2D.Double(i.getX().doubleValue(),i.getY().doubleValue()))
                .toList()
        );
    }

    private void calculate() {
        final int degreeTofit=3;
        final int totalSegmentCountToFit=4;
        System.out.println(String.format("calculating best polyfit of degree %d using %d segments",degreeTofit,totalSegmentCountToFit));

        double globalStartTime = sampleSeries.getMinX();
        double globalEndTime = sampleSeries.getMaxX();
        double avgDelta= (sampleSeries.getMaxX()-sampleSeries.getMinX())/(totalSegmentCountToFit-1);

        List<PolynomialFunction> bestFitResult=List.of();
        List<SegmentSampleData> bestSegments=List.of();
        double bestFitScore = 10000.0;
        for (int k=0;k<10000;k++){
            List<Double> intersectionTimes = new ArrayList<>();
            intersectionTimes.add(globalStartTime);

            double prev = globalStartTime;
            for (int i=0; i<totalSegmentCountToFit-1;i++) {
                double curr = (Math.random() * avgDelta)+prev;
                intersectionTimes.add(curr);
                prev = curr;
            }
            intersectionTimes.add(globalEndTime);

            List<SegmentSampleData> segments = new ArrayList<>();
            for (int i=0; i<totalSegmentCountToFit;i++) {
                var p = getPartition(intersectionTimes.get(i),intersectionTimes.get(i+1),sampleSeries);
                SegmentSampleData segment = getSegmentSampleData(p);
                segments.add(segment);
            }
            PolyfitDto dto = new PolyfitDto(segments, 3);
            CustomPolyFitter customPolyFitter = new CustomPolyFitter(dto);
            List<List<Double>> coeffs = customPolyFitter.calculateOptimalCoeffs();

            List<PolynomialFunction> result = new ArrayList<>();
            for (int i=0; i<totalSegmentCountToFit;i++) {
                PolynomialFunction polynomial = new PolynomialFunction(coeffs.get(i).stream().mapToDouble(Double::doubleValue).toArray());
                result.add(polynomial);
            }
            if (CustomPolyFitter.best<bestFitScore) {
                bestFitScore=CustomPolyFitter.best;
                System.out.println(String.format("New best segmentation at: %s , fit-value= %f", intersectionTimes,bestFitScore));
                bestSegments=segments;
                bestFitResult=result;
            }
        }


        List<Color> colors = List.of(Color.MAGENTA, Color.ORANGE, Color.RED, Color.BLUE, Color.GREEN);

        for (int i=0; i<totalSegmentCountToFit;i++) {
            plotPolynomial(bestSegments.get(0).samples().stream().mapToDouble(s->s.getX()).toArray(), bestFitResult.get(i), plot, 2+i, "p"+String.valueOf(i), colors.get(i), bestSegments.get(i).startTime(), bestSegments.get(i).endTime());
        }

    }

    private SegmentSampleData getSegmentSampleData(PartitionedSample s) {
        List<WeightedObservedPoint> samples = new ArrayList<>();
        for (int i = 0; i < s.points().size(); i++) {
            samples.add(new WeightedObservedPoint(1.0, s.points().get(i).getX(), s.points().get(i).getY()));
        }
        return new SegmentSampleData(s.startTime(), s.endTime(), samples, 30, 10);

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


}