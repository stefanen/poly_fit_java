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

    public int continuityWeight = 30;
    public int derivativeContinuityWeight = 10;
    public int degreeTofit = 3;
    public int totalSegmentCountToFit = 4;
    public static final int SEGMENT_INTERVAL_SEARCH_SPACE_SIZE = 3000;
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
        plot.getDomainAxis().setRange(0.0, 10.0);   // X axis
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
        topPanel.add(new JLabel("Number of segments (1-5):"));
        JTextField segmentCountInput = new JTextField(2);
        segmentCountInput.setText("4");
        segmentCountInput.addActionListener(e -> {
            String text = segmentCountInput.getText();
            this.totalSegmentCountToFit=Integer.parseInt(text);
            recalculatePlot();
        });
        topPanel.add(segmentCountInput);

        topPanel.add(new JLabel("Degree to fit (1-9)"));
        JTextField degreeInput = new JTextField(2);
        degreeInput.setText("3");
        degreeInput.addActionListener(e -> {
            String text = degreeInput.getText();
            this.degreeTofit=Integer.parseInt(text);
            recalculatePlot();
        });
        topPanel.add(degreeInput);

        topPanel.add(new JLabel("Continuity Weight (0-N)"));
        JTextField continuityInput = new JTextField(2);
        continuityInput.setText("30");
        continuityInput.addActionListener(e -> {
            String text = degreeInput.getText();
            this.continuityWeight=Integer.parseInt(text);
            recalculatePlot();
        });
        topPanel.add(continuityInput);

        topPanel.add(new JLabel("Derivative Continuity Weight (0-N)"));
        JTextField derivContinuityInput = new JTextField(2);
        derivContinuityInput.setText("10");
        derivContinuityInput.addActionListener(e -> {
            String text = degreeInput.getText();
            this.derivativeContinuityWeight=Integer.parseInt(text);
            recalculatePlot();
        });
        topPanel.add(derivContinuityInput);

        JButton zoomToggle = new JButton("Toggle edit mode");
        zoomToggle.addActionListener(e -> {
            isReadOnlyMode=!isReadOnlyMode;
            chartPanel.setMouseZoomable(isReadOnlyMode);
        });
        topPanel.add(zoomToggle);

        root.add(topPanel, BorderLayout.NORTH);
        root.add(chartPanel, BorderLayout.CENTER);
        setContentPane(root);
        setSize(1200, 800);
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        setLocationRelativeTo(null);
    }

    private void recalculatePlot() {
        for (int i = plot.getDatasetCount() - 1; i >= 1; i--) {
            plot.setDataset(i, null);
        }
        calculate();
    }

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

        Point2D.Double point = new Point2D.Double(x, y);
        drawnPoints.add(point);
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
        double bestFitScore = 10000.0;
        for (int k = 0; k< SEGMENT_INTERVAL_SEARCH_SPACE_SIZE; k++){

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
        List<Color> colors = List.of(Color.MAGENTA, Color.ORANGE, Color.RED, Color.BLUE, Color.GREEN);

        for (int i = 0; i< totalSegmentCountToFit; i++) {
            plotPolynomial(bestSegments.get(0).samples().stream().mapToDouble(s->s.getX()).toArray(), bestFitResult.get(i), plot, 2+i, "p"+String.valueOf(i)+ " : " +  bestFitResult.get(i).toString().replaceAll("([.][0-9]{3})[0-9]* ","$1 "), colors.get(i), bestSegments.get(i).startTime(), bestSegments.get(i).endTime());
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