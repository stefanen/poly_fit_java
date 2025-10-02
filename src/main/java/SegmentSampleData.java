import org.apache.commons.math3.fitting.WeightedObservedPoint;

import java.util.List;

public record SegmentSampleData(double startTime, double endTime, List<WeightedObservedPoint> samples, double continuityWeight, double derivativeContinuityWeight) {
}
