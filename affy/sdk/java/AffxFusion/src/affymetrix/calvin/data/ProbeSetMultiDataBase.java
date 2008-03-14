package affymetrix.calvin.data;

import java.util.ArrayList;
import java.util.List;

import affymetrix.calvin.parameter.ParameterNameValue;

public class ProbeSetMultiDataBase {

	/** Other metrics associated with the result. */
	private List<ParameterNameValue> metrics = new ArrayList<ParameterNameValue>();

	public List<ParameterNameValue> getMetrics() {
		return metrics;
	}

	public void addMetrics(List<ParameterNameValue> m) {
		this.metrics.addAll(m);
	}

	public void addMetric(ParameterNameValue m) {
		this.metrics.add(m);
	}

	public void clearMetrics() {
		metrics.clear();
	}
}
