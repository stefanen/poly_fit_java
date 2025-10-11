import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

class CustomPolyFitterTest {

    @Test
    void testPolyfitWithoutContinuity() {
        CustomPolyFitter customPolyFitter = getPolyFitDtoForTest(0, 0);
        List<List<Double>> coeffs = customPolyFitter.calculateOptimalCoeffs();
        Assertions.assertArrayEquals(new double[]{0,1,0,0},
                coeffs.get(0).stream().mapToDouble(Double::doubleValue).toArray(),
                1e-5);

        Assertions.assertArrayEquals(new double[]{1,1,0,0},
                coeffs.get(1).stream().mapToDouble(Double::doubleValue).toArray(),
                1e-5);

        double diffAtJunction=Math.abs(new PolynomialFunction(coeffs.get(0).stream().mapToDouble(Double::doubleValue).toArray()).value(10)
                -new PolynomialFunction(coeffs.get(1).stream().mapToDouble(Double::doubleValue).toArray()).value(10));

        Assertions.assertEquals(
                1,
                diffAtJunction,
                1e-5
        );
    }

    @Test
    void testPolyfitWithContinuity() {
        CustomPolyFitter customPolyFitter = getPolyFitDtoForTest(10, 10);
        List<List<Double>> coeffs = customPolyFitter.calculateOptimalCoeffs();
        Assertions.assertArrayEquals(new double[]{-0.036,1.102,-0.038,0.003},
                coeffs.get(0).stream().mapToDouble(Double::doubleValue).toArray(),
                1e-2);

        Assertions.assertArrayEquals(new double[]{-16.328,4.569,-0.240,0.005},
                coeffs.get(1).stream().mapToDouble(Double::doubleValue).toArray(),
                1e-2);

        double diffAtJunction=Math.abs(new PolynomialFunction(coeffs.get(0).stream().mapToDouble(Double::doubleValue).toArray()).value(10)
                -new PolynomialFunction(coeffs.get(1).stream().mapToDouble(Double::doubleValue).toArray()).value(10));
        
        Assertions.assertEquals(
                0.0295,
                diffAtJunction,
                1e-2
        );

        double derivDiffAtJunction=Math.abs(new PolynomialFunction(coeffs.get(0).stream().mapToDouble(Double::doubleValue).toArray()).polynomialDerivative().value(10)
                -new PolynomialFunction(coeffs.get(1).stream().mapToDouble(Double::doubleValue).toArray()).polynomialDerivative().value(10));

        Assertions.assertEquals(
                0.0192,
                derivDiffAtJunction,
                1e-2
        );

    }




    /*
        p_1(t)=t for 0<=t<=10
        p_2(t)=1+t for 10<=t<=20
     */
    private static CustomPolyFitter getPolyFitDtoForTest(double continuityWeight, double derivContinuityWeight) {
        double t1 = 0;
        double t2 = 10;
        double t3 = 20;
        List<WeightedObservedPoint> samples1 = createLinearSamples(t1, t2, 1, 0);
        SegmentSampleData segment1 = new SegmentSampleData(t1, t2, samples1, 0, 0);
        List<WeightedObservedPoint> samples2 = createLinearSamples(t2, t3, 1, 1);
        SegmentSampleData segment2 = new SegmentSampleData(t2, t3, samples2, continuityWeight, derivContinuityWeight);
        List<SegmentSampleData> segments = List.of(segment1, segment2);
        PolyfitDto dto = new PolyfitDto(segments, 3);

        CustomPolyFitter customPolyFitter = new CustomPolyFitter(dto);
        return customPolyFitter;
    }

    /*
        creates samples from y=kx+m (with weight 1), x in [startX,endX]
     */
    private static List<WeightedObservedPoint> createLinearSamples(double startX, double endX, double k, double m) {
        List<WeightedObservedPoint> samples1 = new ArrayList<>();
        for (double i = startX; i < endX; i = i + (endX - startX) / 10) {
            samples1.add(new WeightedObservedPoint(1.0, i, i * k + m));
        }
        return samples1;
    }

}