import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.SimpleVectorValueChecker;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

public class PolyFitterUtil {

    public static record SegmentSampleData(double[] xSamples, double[] ySamples, double startTime, double EndTime) {
    }

    ;

    public static List<List<Double>> polyfitMultiple(List<SegmentSampleData> segments, int coeffCount, double continuityWeight, double derivContinuityWeight) {
        int totalNumberOfSampleDataRows = segments.stream().mapToInt(s -> s.xSamples().length).sum();
        int totalNumberOfCoefficients = segments.size() * coeffCount;
        int numberOfJunctions = segments.size() - 1;
        /*
            We are trying to set up an optimization problem Ac=b, where we want to find c that minimizes |Ac-b|
            A: matrix of dimensions m x n, where first {totalNumberOfSampleDataRows} rows correspond to sample points,
                then, we add rows to enforce continuity. If segment1 and segment2 meet at x0, the condition is p_1(x0)=p_2(x0)
            c: vector of dimensions n x 1, containing all segment-samples' polynomial coefficients
            b: vector of dimensions m x 1, containing all segment-samples' ySamples-values
         */
        int m = totalNumberOfSampleDataRows + numberOfJunctions * 2;
        int n = totalNumberOfCoefficients;

        double[] weights = setupWeights(continuityWeight, derivContinuityWeight, m, totalNumberOfSampleDataRows, numberOfJunctions);
        RealMatrix A = new Array2DRowRealMatrix(m, n);
        RealVector b = new ArrayRealVector(m);
        updateMatrixFromSamplePoints(0, segments, coeffCount, weights, n, A, b);
        updateMatrixWithContinuityConditions(segments, coeffCount, numberOfJunctions, weights, totalNumberOfSampleDataRows, b, n, A, false);
        updateMatrixWithContinuityConditions(segments, coeffCount, numberOfJunctions, weights, totalNumberOfSampleDataRows+numberOfJunctions, b, n, A, true);
        LeastSquaresProblem lsp = new LeastSquaresBuilder()
                .start(new double[n]) // zeroes as guess
                .model(p -> A.operate(p), p -> A.getData()) // linear model
                .target(b)
                .lazyEvaluation(false)
                .maxEvaluations(1000)
                .maxIterations(1000)
                .checkerPair(new SimpleVectorValueChecker(1e-12, 1e-12))
                .build();

        LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
        LeastSquaresOptimizer.Optimum optimizeResult = optimizer.optimize(lsp);

        List<Double> allCoeffs = Arrays.stream(optimizeResult.getPoint().toArray()).boxed().toList();

        return IntStream.range(0, allCoeffs.size() / coeffCount)
                .mapToObj(i -> (List<Double>) new ArrayList<>(allCoeffs.subList(coeffCount * i, coeffCount * (i + 1))))
                .toList();
    }

    private static double[] setupWeights(double continuityWeight, double derivContinuityWeight, int m, int totalNumberOfSampleDataRows, int numberOfJunctions) {
        double[] weights = new double[m];
        for (int i = 0; i < m; i++) {
            if (i < totalNumberOfSampleDataRows) {
                weights[i] = 1.0; //could give weights to samples as input if we wanted
            } else if (i < totalNumberOfSampleDataRows + numberOfJunctions) {
                weights[i] = continuityWeight;
            } else {
                weights[i] = derivContinuityWeight;
            }
        }
        return weights;
    }

    private static void updateMatrixWithContinuityConditions(List<SegmentSampleData> segments, int coeffCount, int numberOfJunctions, double[] weights, int rowIndex, RealVector b, int n, RealMatrix A, boolean isDerivative) {
        int columnIndex = 0;
        for (int i = 0; i < numberOfJunctions; i++) {
            double rowWeight = Math.sqrt(weights[rowIndex]);
            b.setEntry(rowIndex, 0);

            double x0 = segments.get(i).EndTime();
            for (int j = 0; j < n; j++) {
                if (j >= columnIndex && j < columnIndex + coeffCount) {
                    int power = j - columnIndex;
                    if (isDerivative) {
                        A.setEntry(rowIndex, j, power * Math.pow(x0, power - 1) * rowWeight);
                    } else {
                        A.setEntry(rowIndex, j, Math.pow(x0, power) * rowWeight);
                    }
                } else if (j >= columnIndex + coeffCount && j < columnIndex + coeffCount * 2) {
                    int power = j - columnIndex - coeffCount;
                    if (isDerivative) {
                        A.setEntry(rowIndex, j, -1.0 * power * Math.pow(x0, power - 1) * rowWeight);
                    } else {
                        A.setEntry(rowIndex, j, -1.0 * Math.pow(x0, power) * rowWeight);
                    }
                } else {
                    A.setEntry(rowIndex, j, 0.0);
                }
            }
            columnIndex = columnIndex + coeffCount;
            rowIndex++;
        }
    }

    private static void updateMatrixFromSamplePoints(int rowIndex, List<SegmentSampleData> segments, int coeffCount, double[] weights, int n, RealMatrix A, RealVector b) {
        int columnIndex = 0;
        for (var segmentData : segments) {
            for (int i = 0; i < segmentData.xSamples().length; i++) {
                double currX = segmentData.xSamples()[i];
                double currY = segmentData.ySamples()[i];
                double rowWeight = Math.sqrt(weights[rowIndex]);

                for (int j = 0; j < n; j++) {
                    if (j >= columnIndex && j < columnIndex + coeffCount) {
                        A.setEntry(rowIndex, j, Math.pow(currX, j - columnIndex) * rowWeight);
                    } else {
                        A.setEntry(rowIndex, j, 0.0);
                    }
                }
                b.setEntry(rowIndex, currY * rowWeight);
                rowIndex++;
            }
            columnIndex = columnIndex + coeffCount;
        }
    }


    /* fit polynomial p to xSamples-ySamples samples, but also to the derivative having a target value at supplied point p'(x0)=y0  */
    public static double[] fit(double[] x, double[] y, int degree, double x0, double y0, double continuityWeight, double derivContinuityWeight) {
        int nData = x.length;
        int nParams = degree + 1;
        int nRows = nData + 1; //adding a derivative match row

        int derivativeRowIndex = nData;
        double[] weights = new double[nRows];

        for (int i = 0; i < nData; i++) {
            if (x0 == x[i]) {
                weights[i] = continuityWeight;
            } else {
                weights[i] = 1.0;
            }
        }
        weights[derivativeRowIndex] = derivContinuityWeight;

        RealMatrix A = new Array2DRowRealMatrix(nRows, nParams);
        RealVector b = new ArrayRealVector(nRows);

        for (int i = 0; i < nData; i++) {
            double wSqrt = Math.sqrt(weights[i]);
            for (int j = 0; j < nParams; j++) {
                A.setEntry(i, j, Math.pow(x[i], j) * wSqrt);
            }
            b.setEntry(i, y[i] * wSqrt);
        }

        // if p(xSamples)=ax^3+bx^2+cx+d then
        // p'(xSamples)=3ax^2 + 2bx + c = [0,1,2*xSamples,3*xSamples^2] * [ d,c,b,a ]^T
        //so the deriv-row in matrix A should be populated with [0,1,2*xSamples,3*xSamples^2, ...]
        for (int j = 0; j < nParams; j++) {
            if (j == 0) {
                A.setEntry(derivativeRowIndex, j, 0.0);
            } else {
                A.setEntry(derivativeRowIndex, j, j * Math.pow(x0, j - 1) * weights[derivativeRowIndex]);
            }
        }

        b.setEntry(nRows - 1, y0 * weights[derivativeRowIndex]);
        LeastSquaresProblem lsp = new LeastSquaresBuilder()
                .start(new double[nParams]) // initial guess = zeros
                .model(p -> A.operate(p), p -> A.getData()) // linear model
                .target(b)
                .lazyEvaluation(false)
                .maxEvaluations(1000)
                .maxIterations(1000)
                .checkerPair(new SimpleVectorValueChecker(1e-12, 1e-12))
                .build();

        LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
        LeastSquaresOptimizer.Optimum result = optimizer.optimize(lsp);
        return result.getPoint().toArray();
    }

    static Pair<PolynomialFunction, PolynomialFunction> fitTwoAdjoiningPolynomialsSmoothly(double[] p1xSamples, double[] p1ySamples, double[] p2xSamples, double[] p2ySamples, int degree, double cWeight, double dWeight) {

        double commonPointX = p1xSamples[p1ySamples.length - 1];
        double targetDerivativeInCommonPoint = 2; //could be rewritten as a condition that p1'(x_0)=p2'(x_0) but for now we both set a target velocity
        double continuityWeight = cWeight;
        double derivContinuityWeight = dWeight;


        PolynomialFunction p1 = new PolynomialFunction(fit(p1xSamples, p1ySamples, degree, commonPointX, targetDerivativeInCommonPoint, continuityWeight, derivContinuityWeight));
        PolynomialFunction p2 = new PolynomialFunction(fit(p2xSamples, p2ySamples, degree, commonPointX, targetDerivativeInCommonPoint, continuityWeight, derivContinuityWeight));
        System.out.println(Arrays.toString(p1.getCoefficients()));
        System.out.println(Arrays.toString(p2.getCoefficients()));
        return new Pair<>(p1, p2);

    }

    static Pair<PolynomialFunction, PolynomialFunction> fitTwo(double[] p1xSamples, double[] p1ySamples, double[] p2xSamples, double[] p2ySamples, int degree, double cWeight, double dWeight, double s1start, double s1end, double s2end) {


        SegmentSampleData s1 = new SegmentSampleData(p1xSamples, p1ySamples, s1start, s1end );
        SegmentSampleData s2 = new SegmentSampleData(p2xSamples, p2ySamples, s1end, s2end);

        //SegmentSampleData s2 = new SegmentSampleData(p2xSamples,p2ySamples,p1xSamples[p1xSamples.length-1]), p2xSamples[p2xSamples.length-1]));

        List<List<Double>> coeffs = polyfitMultiple(List.of(s1, s2), degree + 1, cWeight, dWeight);
        PolynomialFunction p1 = new PolynomialFunction(coeffs.get(0).stream().mapToDouble(Double::doubleValue).toArray());
        PolynomialFunction p2 = new PolynomialFunction(coeffs.get(1).stream().mapToDouble(Double::doubleValue).toArray());

        return new Pair<>(p1, p2);

    }

    static List<PolynomialFunction> fitThree(double[] p1xSamples, double[] p1ySamples, double[] p2xSamples, double[] p2ySamples, double[] p3xSamples, double[] p3ySamples, int degree, double cWeight, double dWeight, double s1start, double s1end, double s2end, double s3end) {


        SegmentSampleData s1 = new SegmentSampleData(p1xSamples, p1ySamples, s1start, s1end );
        SegmentSampleData s2 = new SegmentSampleData(p2xSamples, p2ySamples, s1end, s2end);
        SegmentSampleData s3 = new SegmentSampleData(p3xSamples, p3ySamples, s2end, s3end);

        //SegmentSampleData s2 = new SegmentSampleData(p2xSamples,p2ySamples,p1xSamples[p1xSamples.length-1]), p2xSamples[p2xSamples.length-1]));

        List<List<Double>> coeffs = polyfitMultiple(List.of(s1, s2,s3), degree + 1, cWeight, dWeight);
        PolynomialFunction p1 = new PolynomialFunction(coeffs.get(0).stream().mapToDouble(Double::doubleValue).toArray());
        PolynomialFunction p2 = new PolynomialFunction(coeffs.get(1).stream().mapToDouble(Double::doubleValue).toArray());
        PolynomialFunction p3 = new PolynomialFunction(coeffs.get(2).stream().mapToDouble(Double::doubleValue).toArray());

        return List.of(p1,p2,p3);

    }




    /* default polyfit algorithm, for comparison */
    static PolynomialFunction fitSamplesDefault(double[] xSamples, double[] ySamples, int degree) {
        WeightedObservedPoints obs = new WeightedObservedPoints();
        for (int i = 0; i < xSamples.length; i++) {
            obs.add(xSamples[i], ySamples[i]);
        }
        PolynomialCurveFitter fitter = PolynomialCurveFitter.create(degree);
        double[] coeff = fitter.fit(obs.toList());
        System.out.println(Arrays.toString(coeff));
        PolynomialFunction polynomial = new PolynomialFunction(coeff);
        return polynomial;
    }
}
