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

public class FunctionFitterToDrawnLine extends JFrame {

    private final JTextField searchDepthInput;
    private final JTextField derivContinuityInput;
    private final JTextField continuityInput;
    private final JTextField degreeInput;
    private final JTextField segmentCountInput;
    public int continuityWeight = 30;
    public int derivativeContinuityWeight = 10;
    public int degreeTofit = 3;
    public int totalSegmentCountToFit = 4;
    public int segmentIntervalSearchSpaceSize = 3000;
    private final XYSeries sampleSeries;
    private final List<Point2D.Double> drawnPoints = new ArrayList<>();
    private final XYPlot plot;
    private final ChartPanel chartPanel;
    private final JFreeChart chart;
    private ScheduledFuture<?> calculationFuture;

    private boolean isReadOnlyMode=false;

    public FunctionFitterToDrawnLine() {
        super("Drag to Draw Points");

        sampleSeries = new XYSeries("Drawn Points");

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
        plot.getDomainAxis().setRange(0.0, 10.0);   // X axisl
        plot.getRangeAxis().setRange(-10.0, 10.0);
        chartPanel = new ChartPanel(chart);
        chartPanel.setMouseZoomable(isReadOnlyMode); // disable zoom

        chartPanel.addMouseMotionListener(new MouseAdapter() {
            @Override
            public void mouseDragged(MouseEvent e) {
                addPointFromMouse(e, chartPanel, plot);
            }
        });
        chartPanel.addMouseListener(new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                if (isReadOnlyMode) {
                    return;
                }
                sampleSeries.clear();
                drawnPoints.clear();
                addPointFromMouse(e, chartPanel, plot);
            }
        });

        JPanel root = new JPanel(new BorderLayout());

        JPanel topPanel = new JPanel();
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
        searchDepthInput.setText("3000");
        searchDepthInput.addActionListener(e -> {
            recalculatePlot();
        });
        topPanel.add(searchDepthInput);

        JButton zoomToggle = new JButton("Toggle drawing mode");
        zoomToggle.addActionListener(e -> {
            isReadOnlyMode=!isReadOnlyMode;
            chartPanel.setMouseZoomable(isReadOnlyMode);
        });
        topPanel.add(zoomToggle);

        JButton restartButton = new JButton("Restart");
        topPanel.add(restartButton);
        restartButton.addActionListener(e -> {
            dispose();
            SwingUtilities.invokeLater(() -> {
                new FunctionFitterToDrawnLine().setVisible(true);
            });
        });

        root.add(topPanel, BorderLayout.NORTH);
        root.add(chartPanel, BorderLayout.CENTER);
        setContentPane(root);
        setSize(1200, 800);
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        setLocationRelativeTo(null);
    }

    private void recalculatePlot() {
        this.derivativeContinuityWeight=Integer.parseInt(derivContinuityInput.getText());
        this.continuityWeight=Integer.parseInt(continuityInput.getText());
        this.degreeTofit=Integer.parseInt(degreeInput.getText());
        this.segmentIntervalSearchSpaceSize=Integer.parseInt(searchDepthInput.getText());
        this.totalSegmentCountToFit=Integer.parseInt(segmentCountInput.getText());
        for (int i = plot.getDatasetCount() - 1; i >= 1; i--) {
            plot.setDataset(i, null);
        }
        calculate();
    }


    Point2D.Double lastSmooth = null;
    double alpha = 0.25;
    private void addPointFromMouse(MouseEvent e, ChartPanel chartPanel, XYPlot plot) {
        if (isReadOnlyMode) {
            return;
        }
        Rectangle2D dataArea = chartPanel.getScreenDataArea();

        ValueAxis xAxis = plot.getDomainAxis();
        ValueAxis yAxis = plot.getRangeAxis();

        double x = xAxis.java2DToValue(
                e.getX(), dataArea, plot.getDomainAxisEdge());
        double y = yAxis.java2DToValue(
                e.getY(), dataArea, plot.getRangeAxisEdge());


        Point2D.Double raw = new Point2D.Double(x, y);

        if (lastSmooth == null) {
            lastSmooth = raw;
        } else {
            x = alpha * raw.x + (1 - alpha) * lastSmooth.x;
            y = alpha * raw.y + (1 - alpha) * lastSmooth.y;
            lastSmooth = new Point2D.Double(x, y);
        }

        drawnPoints.add(raw);

        sampleSeries.add(x, y);



        ScheduledExecutorService scheduledExecutorService = Executors.newScheduledThreadPool(1);
        if (calculationFuture != null) {
            calculationFuture.cancel(false);
        }

        calculationFuture =
                scheduledExecutorService.schedule(this::calculate, 200L, TimeUnit.MILLISECONDS);

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
        System.out.println(String.format("calculating best polyfit of degree %d using %d segments", degreeTofit, totalSegmentCountToFit));
        double globalStartTime = sampleSeries.getMinX();
        double globalEndTime = sampleSeries.getMaxX();
        List<PolynomialFunction> bestFitResult=List.of();
        List<SegmentSampleData> bestSegments=List.of();
        double bestFitScore = 1000000.0;
        for (int k = 0; k< segmentIntervalSearchSpaceSize; k++){

            var intersectionTimes = generateRandomIntersectionTimes(globalStartTime,globalEndTime, totalSegmentCountToFit);

            List<SegmentSampleData> segments = new ArrayList<>();
            for (int i = 0; i< totalSegmentCountToFit; i++) {
                var p = getPartition(intersectionTimes.get(i),intersectionTimes.get(i+1),sampleSeries);
                SegmentSampleData segment = getSegmentSampleData(p);
                segments.add(segment);
            }
            PolyfitDto dto = new PolyfitDto(segments, degreeTofit);
            CustomPolyFitter customPolyFitter = new CustomPolyFitter(dto);
            List<List<Double>> coeffs = customPolyFitter.calculateOptimalCoeffs();

            List<PolynomialFunction> result = new ArrayList<>();
            for (int i = 0; i< totalSegmentCountToFit; i++) {
                PolynomialFunction polynomial = new PolynomialFunction(coeffs.get(i).stream().mapToDouble(Double::doubleValue).toArray());
                result.add(polynomial);
            }
            if (CustomPolyFitter.lastRes <bestFitScore) {
                bestFitScore=CustomPolyFitter.lastRes;
                System.out.println(String.format("New best segmentation at: %s , fit-value= %f, after %d tries", intersectionTimes,bestFitScore,k));
                bestSegments=segments;
                bestFitResult=result;
            }
            /*
            if (k%100==0) {
                System.out.println(String.format("segmentation at: %s , fit-value= %f, after %d tries", intersectionTimes,CustomPolyFitter.lastRes,k));

            }*/

        }
        List<Color> colors = List.of(Color.MAGENTA, Color.ORANGE, Color.RED, Color.BLUE, Color.GREEN, Color.WHITE, Color.DARK_GRAY, Color.CYAN, Color.PINK, Color.YELLOW);

        for (int i = 0; i< totalSegmentCountToFit; i++) {
            plotPolynomial(bestSegments.get(0).samples().stream().mapToDouble(s->s.getX()).toArray(), bestFitResult.get(i), plot, 2+i, "p"+String.valueOf(i)+ " : " +  bestFitResult.get(i).toString().replaceAll("([.][0-9]{3})[0-9]*","$1 "), colors.get(i), bestSegments.get(i).startTime(), bestSegments.get(i).endTime());
        }

        chart.setTitle(String.format("Showing best poly-fit to drawn sample-data, using degree=%d and segmentCount=%d", degreeTofit, totalSegmentCountToFit));
    }


    public static List<Double> generateRandomIntersectionTimes(double min, double max, int n) {
        Random rand = new Random();
        List<Double> cuts = new ArrayList<>();
        cuts.add(min);
        for (int i = 1; i < n; i++) {
            cuts.add(min + rand.nextDouble() * (max-min));
        }
        cuts.add(max);


        return cuts.stream().sorted().toList();
    }


    private SegmentSampleData getSegmentSampleData(PartitionedSample s) {
        List<WeightedObservedPoint> samples = new ArrayList<>();
        for (int i = 0; i < s.points().size(); i++) {
            samples.add(new WeightedObservedPoint(1.0, s.points().get(i).getX(), s.points().get(i).getY()));
        }
        return new SegmentSampleData(s.startTime(), s.endTime(), samples, continuityWeight, derivativeContinuityWeight);

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