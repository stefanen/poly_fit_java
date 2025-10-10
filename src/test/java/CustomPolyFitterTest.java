import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.junit.jupiter.params.provider.ValueSources;

import java.util.ArrayList;
import java.util.List;

class CustomPolyFitterTest {


    @ParameterizedTest
    //@ValueSource(ints = {0,1,10})
    @CsvSource({"0,0, 0,1,0,0, 1,1,0,0"
    ,
            "10,10, 0,1,0,0, 1,1,0,0"
    })

    void testPolyFit(int a, int b, int p_1_a0, int p_1_a1, int p_1_a2, int p_1_a3, int p_2_a0,int  p_2_a1, int p_2_a2, int p_2_a3) {
        double t1=0;
        double t2=10;
        double t3=20;
        List<WeightedObservedPoint> samples1 = createLinearSamples(1,1);
        SegmentSampleData segment1 = new SegmentSampleData(t1, t2, samples1, 0,0);
        List<WeightedObservedPoint> samples2 = createLinearSamples(11,12);
        SegmentSampleData segment2 = new SegmentSampleData(t2, t3, samples2, a,b);
        List<SegmentSampleData> segments = List.of(segment1,segment2);
        PolyfitDto dto = new PolyfitDto(segments,3);
        CustomPolyFitter customPolyFitter = new CustomPolyFitter(dto);

        List<List<Double>> coeffs = customPolyFitter.calculateOptimalCoeffs();

        double delta =0;
        //TODO cleanup
        delta+=Math.abs(coeffs.get(0).get(0)-p_1_a0);
        delta+=Math.abs(coeffs.get(0).get(1)-p_1_a1);
        delta+=Math.abs(coeffs.get(0).get(2)-p_1_a2);
        delta+=Math.abs(coeffs.get(0).get(3)-p_1_a3);

        delta+=Math.abs(coeffs.get(1).get(0)-p_2_a0);
        delta+=Math.abs(coeffs.get(1).get(1)-p_2_a1);
        delta+=Math.abs(coeffs.get(1).get(2)-p_2_a2);
        delta+=Math.abs(coeffs.get(1).get(3)-p_2_a3);

        System.out.println(coeffs);

        System.out.println(delta);
        Assertions.assertTrue(delta<0.001);
    }

    private static List<WeightedObservedPoint> createLinearSamples(int startX,int startY) {
        List<WeightedObservedPoint> samples1 = new ArrayList<>();
        for (int i=0; i<4; i++) {
            samples1.add(new WeightedObservedPoint(1.0, startX+i, startY+i));
        }
        return samples1;
    }


}