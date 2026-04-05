package com.stefanen;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleEdge;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.ScheduledFuture;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

public class FunctionFitterToBezierCurve extends JFrame {

    public static final int CONTROL_DATA_SET_INDEX = 10;
    private JTextField searchDepthInput;
    private JTextField derivContinuityInput;
    private JTextField continuityInput;
    private JTextField degreeInput;
    private JTextField segmentCountInput;
    private XYSeries controlDataSeries;
    public int continuityWeight = 30;
    public int derivativeContinuityWeight = 10;
    public int degreeTofit = 3;
    public int totalSegmentCountToFit = 4;
    public int segmentIntervalSearchSpaceSize = 200;

    private int selectedBezierPointIndex = -1;

    private final XYSeries sampleSeries;
    private final XYPlot plot;
    private final ChartPanel chartPanel;
    private final JFreeChart chart;
    private ScheduledFuture<?> calculationFuture;
    private boolean isReadOnlyMode = false;

    public FunctionFitterToBezierCurve() {
        super("Bezier-curve to polyfunction");

        sampleSeries = new XYSeries("Points to fit (from Bezier curve)");

        XYSeriesCollection dataset = new XYSeriesCollection(sampleSeries);

        chart = ChartFactory.createXYLineChart(
                "Polyfit where segment-interval-times are variables",
                "X",
                "Y",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );
        LegendTitle legend = chart.getLegend();
        legend.setPosition(RectangleEdge.RIGHT);
        legend.setItemFont(new Font("SansSerif", Font.PLAIN, 12));


        plot = chart.getXYPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(true, false);
        plot.setRenderer(renderer);
        plot.getDomainAxis().setRange(0.0, 20.0);
        plot.getRangeAxis().setRange(-20.0, 20.0);
        chartPanel = new ChartPanel(chart);
        chartPanel.setMouseZoomable(isReadOnlyMode); // disable zoom

        initBezierControlPoints();

        JPanel root = new JPanel(new BorderLayout());
        JPanel topPanel = new JPanel();

        setupUserControls(topPanel);

        root.add(topPanel, BorderLayout.NORTH);
        root.add(chartPanel, BorderLayout.CENTER);
        setContentPane(root);
        setSize(1200, 800);
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        setLocationRelativeTo(null);
    }

    private void setupUserControls(JPanel topPanel) {
        topPanel.add(new JLabel("#Segments to fit (1-9):"));
        segmentCountInput = new JTextField(2);
        segmentCountInput.setText("4");
        segmentCountInput.addActionListener(e -> {
            recalculatePlot();
        });
        topPanel.add(segmentCountInput);

        topPanel.add(new JLabel("Degree to fit"));
        degreeInput = new JTextField(2);
        degreeInput.setText("3");
        degreeInput.addActionListener(e -> {
            recalculatePlot();
        });
        topPanel.add(degreeInput);

        topPanel.add(new JLabel("Continuity Weight"));
        continuityInput = new JTextField(4);
        continuityInput.setText("30");
        continuityInput.addActionListener(e -> {
            recalculatePlot();
        });
        topPanel.add(continuityInput);

        topPanel.add(new JLabel("Derivative Continuity Weight"));
        derivContinuityInput = new JTextField(4);
        derivContinuityInput.setText("10");
        derivContinuityInput.addActionListener(e -> {
            recalculatePlot();
        });
        topPanel.add(derivContinuityInput);

        topPanel.add(new JLabel("Fitting search depth"));
        searchDepthInput = new JTextField(7);
        searchDepthInput.setText("200");
        searchDepthInput.addActionListener(e -> {
            recalculatePlot();
        });
        topPanel.add(searchDepthInput);

        JButton zoomToggle = new JButton("Toggle drawing mode");
        zoomToggle.addActionListener(e -> {
            isReadOnlyMode = !isReadOnlyMode;
            chartPanel.setMouseZoomable(isReadOnlyMode);
        });
        topPanel.add(zoomToggle);

        JButton restartButton = new JButton("Restart");
        topPanel.add(restartButton);
        restartButton.addActionListener(e -> {
            dispose();
            SwingUtilities.invokeLater(() -> {
                new FunctionFitterToBezierCurve().setVisible(true);
            });
        });
    }

    private void setupControlPointsMouseHandlers() {
        chartPanel.addMouseMotionListener(new MouseAdapter() {
            @Override
            public void mouseDragged(MouseEvent e) {
                if (selectedBezierPointIndex == -1) {
                    return;
                }
                Rectangle2D dataArea = chartPanel.getScreenDataArea();

                ValueAxis xAxis = plot.getDomainAxis();
                ValueAxis yAxis = plot.getRangeAxis();

                double x = xAxis.java2DToValue(
                        e.getX(), dataArea, plot.getDomainAxisEdge());
                double y = yAxis.java2DToValue(
                        e.getY(), dataArea, plot.getRangeAxisEdge());


                List<XYDataItem> controlItems = controlDataSeries.getItems().stream().map(i -> (XYDataItem) i).toList();
                var allControlPoints = controlItems.stream().map((XYDataItem i) -> new Point2D.Double(i.getX().doubleValue(), i.getY().doubleValue())).toList();

                controlDataSeries.clear();
                int i = 0;
                for (var controlPoint : allControlPoints) {
                    if (i == selectedBezierPointIndex) {
                        controlDataSeries.add(x, y);
                    } else {
                        controlDataSeries.add(controlPoint.x, controlPoint.y);
                    }
                    i++;
                }


                calculateDelayed();
            }
        });

        chartPanel.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseReleased(MouseEvent e) {
                selectedBezierPointIndex = -1;
            }
        });

        chartPanel.addMouseListener(new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                if (isReadOnlyMode) {
                    return;
                }
                sampleSeries.clear();
                Rectangle2D dataArea = chartPanel.getScreenDataArea();

                ValueAxis xAxis = plot.getDomainAxis();
                ValueAxis yAxis = plot.getRangeAxis();

                double x = xAxis.java2DToValue(
                        e.getX(), dataArea, plot.getDomainAxisEdge());
                double y = yAxis.java2DToValue(
                        e.getY(), dataArea, plot.getRangeAxisEdge());

                double bestDistance = Integer.MAX_VALUE; // depends on axis scale
                List<XYDataItem> controlItems = controlDataSeries.getItems().stream().map(i -> (XYDataItem) i).toList();
                var controlPoints = controlItems.stream().map((XYDataItem i) -> new Point2D.Double(i.getX().doubleValue(), i.getY().doubleValue())).toList();

                int i = 0;
                for (Point2D.Double p : controlPoints) {
                    double dx = p.x - x;
                    double dy = p.y - y;

                    double curr = Math.sqrt(dx * dx + dy * dy);
                    if (curr < bestDistance) {
                        bestDistance = curr;
                        selectedBezierPointIndex = i;
                    }
                    i++;
                }

            }
        });
    }

    private void initBezierControlPoints() {
        final XYSeries control;
        Point2D.Double p0 = new Point2D.Double(0, 0);
        Point2D.Double p1 = new Point2D.Double(2, 5);
        Point2D.Double p2 = new Point2D.Double(5, 5);

        Point2D.Double p3 = new Point2D.Double(8, 0);

        Point2D.Double p4 = new Point2D.Double(11, -5);
        Point2D.Double p5 = new Point2D.Double(12, 10);
        Point2D.Double p6 = new Point2D.Double(16, 5);

        control = new XYSeries("Control points", false);
        control.add(p0.x, p0.y);
        control.add(p1.x, p1.y);
        control.add(p2.x, p2.y);
        control.add(p3.x, p3.y);
        control.add(p4.x, p4.y);
        control.add(p5.x, p5.y);
        control.add(p6.x, p6.y);
        this.controlDataSeries = control;

        XYSeriesCollection datasetControl = new XYSeriesCollection(controlDataSeries);
        XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer(true, true);
        renderer2.setSeriesStroke(0,
                new BasicStroke(
                        2.0f,
                        BasicStroke.CAP_ROUND,
                        BasicStroke.JOIN_ROUND,
                        1.0f,
                        new float[]{6.0f, 4.0f},
                        0.0f));
        renderer2.setSeriesPaint(0, Color.GRAY);
        plot.setDataset(CONTROL_DATA_SET_INDEX, datasetControl);
        plot.setRenderer(CONTROL_DATA_SET_INDEX, renderer2);
        plot.setDatasetRenderingOrder(DatasetRenderingOrder.REVERSE);
        updateSampleSeriesFromControlPoints();
        setupControlPointsMouseHandlers();
        calculateDelayed();
    }

    public record BestResult(List<Double> intersectionTimes, List<SegmentSampleData> segments,
                             List<PolynomialFunction> polynomials, double score) {

    }

    private void calculate() {
        updateSampleSeriesFromControlPoints();
        System.out.println(String.format("calculating best polyfit of degree %d using %d segments", degreeTofit, totalSegmentCountToFit));
        double globalStartTime = sampleSeries.getMinX();
        double globalEndTime = sampleSeries.getMaxX();

        BestResult bestResult = new BestResult(List.of(), List.of(), List.of(), 1000000.0);
        bestResult = getBestSegmentationFromRandomCuts(globalStartTime, globalEndTime, bestResult);
        bestResult = optimizeBestResultWithSmallPertubations(bestResult);

        List<Color> colors = List.of(Color.MAGENTA, Color.ORANGE, Color.RED, Color.BLUE, Color.GREEN, Color.WHITE, Color.DARK_GRAY, Color.CYAN, Color.PINK, Color.YELLOW);

        for (int i = 0; i < totalSegmentCountToFit; i++) {
            plotPolynomial(bestResult.segments().get(0).samples().stream().mapToDouble(s -> s.getX()).toArray(), bestResult.polynomials().get(i), plot, 2 + i, "p" + String.valueOf(i) + " : " + bestResult.polynomials().get(i).toString().replaceAll("([.][0-9]{3})[0-9]*", "$1 "), colors.get(i), bestResult.segments().get(i).startTime(), bestResult.segments().get(i).endTime());
        }

        double startY = 0;
        double endY = 0;
        double startDY_DX = 0;
        double endDY_DX = 0;
        for (int j = 0; j < bestResult.polynomials().size(); j++) {
            startY = bestResult.polynomials().get(j).value(bestResult.intersectionTimes().get(j));
            startDY_DX = bestResult.polynomials().get(j).polynomialDerivative().value(bestResult.intersectionTimes().get(j));
            if (j > 0) {
                System.out.println(String.format("For Segment %d: y diff = %f, y' diff = %f", j, endY - startY, endDY_DX - startDY_DX));
            }
            endY = bestResult.polynomials().get(j).value(bestResult.intersectionTimes().get(j + 1));
            endDY_DX = bestResult.polynomials().get(j).polynomialDerivative().value(bestResult.intersectionTimes().get(j + 1));
            //System.out.println(String.format("%f %f %f %f",startY,startDY_DX,endY,endDY_DX));

        }

        chart.setTitle(String.format("Showing best segmentation and best polynomial-fit to drawn Bezier-curve. \nUsing degree=%d and segmentCount=%d. Fit-score=%f", degreeTofit, totalSegmentCountToFit, bestResult.score()));
    }

    private BestResult optimizeBestResultWithSmallPertubations(BestResult bestResult) {
        List<Double> intersectionTimes = bestResult.intersectionTimes();
        boolean goLeft = true;
        for (int j = 1; j < intersectionTimes.size() - 1; j++) {
            for (int k = 0; k < 1000; k++) {
                var old = intersectionTimes.get(j);
                double step = 0.01;
                if (goLeft) {
                    if (intersectionTimes.get(j - 1) < (old - step)) {
                        intersectionTimes.set(j, old - step);
                    }
                } else {
                    if (intersectionTimes.get(j + 1) > (old + step)) {
                        intersectionTimes.set(j, old + step);
                    }
                }
                List<SegmentSampleData> segments = new ArrayList<>();
                for (int i = 0; i < totalSegmentCountToFit; i++) {
                    var p = getPartition(intersectionTimes.get(i), intersectionTimes.get(i + 1), sampleSeries);
                    SegmentSampleData segment = getSegmentSampleData(p);
                    segments.add(segment);
                }
                PolyfitDto dto = new PolyfitDto(segments, degreeTofit);
                CustomPolyFitter customPolyFitter = new CustomPolyFitter(dto);
                List<List<Double>> coeffs = customPolyFitter.calculateOptimalCoeffs();

                List<PolynomialFunction> result = new ArrayList<>();
                for (int i = 0; i < totalSegmentCountToFit; i++) {
                    PolynomialFunction polynomial = new PolynomialFunction(coeffs.get(i).stream().mapToDouble(Double::doubleValue).toArray());
                    result.add(polynomial);
                }
                //System.out.println("Trying:"+intersectionTimes + ":"+CustomPolyFitter.lastRes);
                if (CustomPolyFitter.lastRes < bestResult.score()) {
                    bestResult = new BestResult(intersectionTimes, segments, result, CustomPolyFitter.lastRes);
                    System.out.println(String.format("New best resegmentation at: %s , fit-value= %f, after %d tries", intersectionTimes, bestResult.score(), k));
                } else {
                    if (k > 2) {
                        break;
                    }
                    intersectionTimes.set(j, old);
                    goLeft = !goLeft;
                }

            }
        }
        return bestResult;
    }

    private BestResult getBestSegmentationFromRandomCuts(double globalStartTime, double globalEndTime, BestResult bestResult) {
        List<Double> intersectionTimes;
        //get a sensible starting point for intersectionTimes by random cuts
        for (int k = 0; k < segmentIntervalSearchSpaceSize; k++) {

            intersectionTimes = generateRandomIntersectionTimes(globalStartTime, globalEndTime, totalSegmentCountToFit);

            List<SegmentSampleData> segments = new ArrayList<>();
            for (int i = 0; i < totalSegmentCountToFit; i++) {
                var p = getPartition(intersectionTimes.get(i), intersectionTimes.get(i + 1), sampleSeries);
                SegmentSampleData segment = getSegmentSampleData(p);
                segments.add(segment);
            }
            PolyfitDto dto = new PolyfitDto(segments, degreeTofit);
            CustomPolyFitter customPolyFitter = new CustomPolyFitter(dto);
            List<List<Double>> coeffs = customPolyFitter.calculateOptimalCoeffs();

            List<PolynomialFunction> result = new ArrayList<>();
            for (int i = 0; i < totalSegmentCountToFit; i++) {
                PolynomialFunction polynomial = new PolynomialFunction(coeffs.get(i).stream().mapToDouble(Double::doubleValue).toArray());
                result.add(polynomial);
            }
            if (CustomPolyFitter.lastRes < bestResult.score()) {
                bestResult = new BestResult(intersectionTimes, segments, result, CustomPolyFitter.lastRes);
                System.out.println(String.format("New best segmentation at: %s , fit-value= %f, after %d tries", intersectionTimes, bestResult.score(), k));
            }
        }
        return bestResult;
    }

    private void calculateDelayed() {
        ScheduledExecutorService scheduledExecutorService = Executors.newScheduledThreadPool(1);
        if (calculationFuture != null) {
            calculationFuture.cancel(false);
        }

        calculationFuture =
                scheduledExecutorService.schedule(this::calculate, 200L, TimeUnit.MILLISECONDS);
    }

    private void recalculatePlot() {
        this.derivativeContinuityWeight = Integer.parseInt(derivContinuityInput.getText());
        this.continuityWeight = Integer.parseInt(continuityInput.getText());
        this.degreeTofit = Integer.parseInt(degreeInput.getText());
        this.segmentIntervalSearchSpaceSize = Integer.parseInt(searchDepthInput.getText());
        this.totalSegmentCountToFit = Integer.parseInt(segmentCountInput.getText());
        for (int i = plot.getDatasetCount() - 1; i >= 1; i--) {
            if (i != CONTROL_DATA_SET_INDEX) {
                plot.setDataset(i, null);
            }
        }
        calculate();
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            new FunctionFitterToBezierCurve().setVisible(true);
        });
    }

    record PartitionedSample(double startTime, double endTime, List<Point2D.Double> points) {

    }

    private PartitionedSample getPartition(double startTime, double endTime, XYSeries allSamples) {
        List<XYDataItem> items = sampleSeries.getItems().stream().map(i -> (XYDataItem) i).toList();

        return new PartitionedSample(startTime, endTime, items.stream()
                .filter((XYDataItem i) -> i.getX().doubleValue() >= startTime && i.getX().doubleValue() < endTime)
                .map((XYDataItem i) -> new Point2D.Double(i.getX().doubleValue(), i.getY().doubleValue()))
                .toList()
        );
    }

    private void updateSampleSeriesFromControlPoints() {
        sampleSeries.clear();

        List<XYDataItem> controlItems = controlDataSeries.getItems().stream().map(i -> (XYDataItem) i).toList();
        var controlPoints = controlItems.stream().map((XYDataItem i) -> new Point2D.Double(i.getX().doubleValue(), i.getY().doubleValue())).toList();

        calculateSamplesOnBezierCurve(controlPoints, 0);
        calculateSamplesOnBezierCurve(controlPoints, 3);


    }

    private void calculateSamplesOnBezierCurve(List<Point2D.Double> controlPoints, int startIndex) {
        var p0 = controlPoints.get(startIndex);
        var p1 = controlPoints.get(startIndex + 1);
        var p2 = controlPoints.get(startIndex + 2);
        var p3 = controlPoints.get(startIndex + 3);

        for (int i = 0; i <= 200; i++) {
            double t = i / 200.0;

            double u = 1 - t;

            double x =
                    u * u * u * p0.x +
                            3 * u * u * t * p1.x +
                            3 * u * t * t * p2.x +
                            t * t * t * p3.x;
            double y =
                    u * u * u * p0.y +
                            3 * u * u * t * p1.y +
                            3 * u * t * t * p2.y +
                            t * t * t * p3.y;

            sampleSeries.add(x, y);
        }
    }


    public static List<Double> generateRandomIntersectionTimes(double min, double max, int n) {
        Random rand = new Random();
        List<Double> cuts = new ArrayList<>();
        cuts.add(min);
        for (int i = 1; i < n; i++) {
            cuts.add(min + rand.nextDouble() * (max - min));
        }
        cuts.add(max);


        return cuts.stream().sorted().collect(Collectors.toList());
    }


    private SegmentSampleData getSegmentSampleData(PartitionedSample s) {
        List<WeightedObservedPoint> samples = new ArrayList<>();
        for (int i = 0; i < s.points().size(); i++) {
            samples.add(new WeightedObservedPoint(1.0, s.points().get(i).getX(), s.points().get(i).getY()));
        }
        return new SegmentSampleData(s.startTime(), s.endTime(), samples, continuityWeight, derivativeContinuityWeight);

    }


    private static void plotPolynomial(double[] x, PolynomialFunction polynomial, XYPlot plot, int id, String name, Color color, double start, double end) {
        List<Double> xValues = new ArrayList<>();
        List<Double> yValues = new ArrayList<>();

        XYSeries series2 = new XYSeries(name);
        for (double i = start; i < end; i += 0.01) {

            series2.add(i, polynomial.value(i));
            xValues.add(i);
            yValues.add(polynomial.value(i));
        }
        XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer(true, true);
        renderer2.setSeriesPaint(0, color);
        plot.setDataset(id, new XYSeriesCollection(series2));
        plot.setRenderer(id, renderer2);
        plot.setDatasetRenderingOrder(DatasetRenderingOrder.REVERSE);
    }


}