package com.stefanen;

import java.util.List;

public record PolyfitDto(List<SegmentSampleData> segments, int degreeToFit) {
}
