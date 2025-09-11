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
import org.apache.commons.math3.optim.PointVectorValuePair;
import org.apache.commons.math3.optim.SimpleVectorValueChecker;
import org.apache.commons.math3.util.Pair;

import java.util.Arrays;

public class PolyFitterUtil {





    /* fit polynomial p to x-y samples, but also to the derivative having a target value at supplied point p'(x0)=y0  */
    public static double[] fit(double[] x, double[] y, int degree, double x0, double y0, double continuityWeight, double derivContinuityWeight) {
        int nData = x.length;
        int nParams = degree + 1;
        int nRows = nData +1; //adding a derivative match row

        int derivativeRowIndex=nData;
        double[] weights = new double[nRows];

        for (int i = 0; i < nData; i++) {
            if (x0 == x[i]) {
                weights[i] = continuityWeight;
            } else {
                weights[i] = 1.0;
            }
        }
        weights[derivativeRowIndex]=derivContinuityWeight;

        RealMatrix A = new Array2DRowRealMatrix(nRows, nParams);
        RealVector b = new ArrayRealVector(nRows);

        for (int i = 0; i < nData; i++) {
            double wSqrt = Math.sqrt(weights[i]);
            for (int j = 0; j < nParams; j++) {
                A.setEntry(i, j, Math.pow(x[i], j)*wSqrt);
            }
            b.setEntry(i, y[i]*wSqrt);
        }

        for (int j = 0; j < nParams; j++) {
            if (j == 0) {
                A.setEntry(derivativeRowIndex, j, 0.0);
            } else {
                A.setEntry(derivativeRowIndex, j, j * Math.pow(x0, j - 1)*weights[derivativeRowIndex]);
            }
        }

        b.setEntry(nRows-1, y0*weights[derivativeRowIndex]);

        // Least squares problem
        LeastSquaresProblem lsp = new LeastSquaresBuilder()
                .start(new double[nParams]) // initial guess = zeros
                .model(p -> A.operate(p), p -> A.getData()) // linear model
                .target(b)
                .lazyEvaluation(false)
                .maxEvaluations(1000)
                .maxIterations(1000)
                .checkerPair(new SimpleVectorValueChecker(1e-12, 1e-12))
                .build();

        // Solve with Levenberg–Marquardt (works fine for linear systems too)
        LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
        LeastSquaresOptimizer.Optimum result = optimizer.optimize(lsp);
        return result.getPoint().toArray();
    }



    static Pair<PolynomialFunction, PolynomialFunction> fitTwoAdjoiningPolynomialsSmoothly(double[] p1xSamples, double[] p1ySamples, double[] p2xSamples, double[] p2ySamples, int degree) {

        double commonPointX=p1xSamples[p1ySamples.length-1];
        double targetDerivativeInCommonPoint = 2; //could be rewritten as a condition that p1'(x_0)=p2'(x_0) but for now we both set a target velocity
        double continuityWeight = 10;
        double derivContinuityWeight = 2;

        PolynomialFunction p1 = new PolynomialFunction(fit(p1xSamples,p1ySamples,degree,commonPointX,targetDerivativeInCommonPoint,continuityWeight,derivContinuityWeight));
        PolynomialFunction p2 = new PolynomialFunction(fit(p2xSamples,p2ySamples,degree,commonPointX,targetDerivativeInCommonPoint,continuityWeight, derivContinuityWeight));
        System.out.println(Arrays.toString( p1.getCoefficients()));
        System.out.println(Arrays.toString( p2.getCoefficients()));
        return new Pair<>(p1,p2);

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
