import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.SimpleVectorValueChecker;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

public class CustomPolyFitter {
    private final PolyfitDto polyfitDto;
    private final int totalNumberOfSampleDataRows;
    private final int numberOfJunctions;
    private final int m;
    private final int n;
    private final int coeffCount;
    private final RealMatrix A;
    private final RealVector b;
    private final double[] weights;

    public CustomPolyFitter(PolyfitDto polyfitDto) {
        this.polyfitDto = polyfitDto;
        totalNumberOfSampleDataRows = polyfitDto.segments().stream().mapToInt(s -> s.samples().size()).sum();
        coeffCount = polyfitDto.degreeToFit() + 1;
        n = polyfitDto.segments().size() * coeffCount;
        numberOfJunctions = polyfitDto.segments().size() - 1;
        m = totalNumberOfSampleDataRows + numberOfJunctions * 2;
        /*
        In order to find best fitting polynomial coefficients (for all segments) we translate it into an optimization problem
        of finding c that minimizes |Ac-b| (c being the vector of all coefficients)
        */
        A = new Array2DRowRealMatrix(m, n);
        b = new ArrayRealVector(m);
        weights = setupWeights();
        addConditionsFromSamplePoints();
        addContinuityConditions();
    }

    /*
     * In the returned list, list.get(i) is a list of size {coeffCount} containing the optimal polynomial coefficients for the i:th segment passed in
     *
     * Example: Suppose we have 3 segments defined by intervals [t_0,t_1], [t_1,t_2], [t_2,t_3], and a handful of sample points {(x_i,y_i)} in each segment,
     * this algorithm then tries to find 3 polynomials p_1(t), p_2(t), p_3(t) that best matches the sample points, AND that best matches the continuity conditions
     *  p_1(t_1)=p_2(t_1)   (continuity from segment 1 to 2)
     *  p_2(t_2)=p_3(t_2)   (continuity from segment 2 to 3)
     *  p_1'(t_1)=p_2'(t_1) (derivative-continuity from segment 1 to 2)
     *  p_2'(t_2)=p_3'(t_2) (derivative-continuity from segment 2 to 3)
     *
     */
    public List<List<Double>> calculateOptimalCoeffs() {
        LeastSquaresProblem lsp = new LeastSquaresBuilder()
                .start(new double[n])
                .model(p -> A.operate(p), p -> A.getData())
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

    private double[] setupWeights() {
        double[] weights = new double[m];
        int i = 0;
        for (var segment : polyfitDto.segments()) {
            for (var sample : segment.samples()) {
                weights[i] = sample.getWeight();
                i++;
            }
        }
        for (int j = 0; j < numberOfJunctions; j++) {
            weights[i + j] = polyfitDto.segments().get(j + 1).continuityWeight();
            weights[i + numberOfJunctions + j] = polyfitDto.segments().get(j + 1).derivativeContinuityWeight();
        }
        return weights;
    }

    private void addContinuityConditions() {
        int rowIndex = totalNumberOfSampleDataRows;
        int columnIndex = 0;
        for (int i = 0; i < numberOfJunctions; i++) {
            double rowWeight = Math.sqrt(weights[rowIndex]);
            double rowWeightDerivative = Math.sqrt(weights[rowIndex+numberOfJunctions]);
            b.setEntry(rowIndex, 0);
            b.setEntry(rowIndex+numberOfJunctions, 0);
            double x0 = polyfitDto.segments().get(i).endTime();
            for (int j = 0; j < n; j++) {
                if (j >= columnIndex && j < columnIndex + coeffCount) {
                    int power = j - columnIndex;
                    A.setEntry(rowIndex+numberOfJunctions, j, power * Math.pow(x0, power - 1) * rowWeightDerivative);
                    A.setEntry(rowIndex, j, Math.pow(x0, power) * rowWeight);
                } else if (j >= columnIndex + coeffCount && j < columnIndex + coeffCount * 2) {
                    int power = j - columnIndex - coeffCount;
                    A.setEntry(rowIndex+numberOfJunctions, j, -1.0 * power * Math.pow(x0, power - 1) * rowWeightDerivative);
                    A.setEntry(rowIndex, j, -1.0 * Math.pow(x0, power) * rowWeight);
                } else {
                    A.setEntry(rowIndex, j, 0.0);
                    A.setEntry(rowIndex+numberOfJunctions, j, 0.0);
                }
            }
            columnIndex = columnIndex + coeffCount;
            rowIndex++;
        }
    }

    private void addConditionsFromSamplePoints() {
        int rowIndex = 0;
        int columnIndex = 0;
        for (var segmentData : polyfitDto.segments()) {
            for (int i = 0; i < segmentData.samples().size(); i++) {
                double currX = segmentData.samples().get(i).getX();
                double currY = segmentData.samples().get(i).getY();
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

}
