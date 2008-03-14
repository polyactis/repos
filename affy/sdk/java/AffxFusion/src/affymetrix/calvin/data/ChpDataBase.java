package affymetrix.calvin.data;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import affymetrix.calvin.parameter.ParameterNameValue;

public class ChpDataBase {

	/** The generic data item. */
	protected GenericData genericData = new GenericData();

	/**
	 * gets the parameters
	 * 
	 * @return The parameters.
	 */
	protected List<ParameterNameValue> getParams(String type) {
		List<ParameterNameValue> algParams = new ArrayList<ParameterNameValue>();
		List<ParameterNameValue> allParams = genericData.getHeader().getGenericDataHdr().getNameValParams();
		for (int i = 0; i < allParams.size(); i++) {
			ParameterNameValue param = allParams.get(i);
			String name = param.getName();
			if (name.startsWith(type)) {
				ParameterNameValue algParam = new ParameterNameValue(param);
				algParam.setName(name.substring(type.length(), name.length()));
				algParams.add(algParam);
			}
		}
		return algParams;
	}

	/**
	 * gets a parameter value as a string.
	 * 
	 * @param name
	 *          The name of the parameter.
	 * @return The string representation.
	 */
	protected String getStringFromGenericHdr(String name) {
		GenericDataHeader hdr = genericData.getHeader().getGenericDataHdr();
		if (hdr == null) {
			return null;
		}
		ParameterNameValue paramType = hdr.findNameValParam(name);
		if (paramType != null) {
			return paramType.getValueText();
		}
		return null;
	}

	protected void setStringToGenericHdr(String name, String value) {
		setStringToGenericHdr(name, value, 0);
	}

	protected void setStringToGenericHdr(String name, String value, int reserve) {
		ParameterNameValue paramType = new ParameterNameValue();
		paramType.setName(name);
		paramType.setValueText(value, reserve);
		GenericDataHeader hdr = genericData.getHeader().getGenericDataHdr();
		hdr.addNameValParam(paramType);
	}

	/**
	 * Gets a parameter value as an integer.
	 * 
	 * @param name
	 *          The name of the parameter.
	 * @return The integer representation.
	 */
	protected int getIntFromGenericHdr(String name) {
		int result = 0;
		GenericDataHeader hdr = genericData.getHeader().getGenericDataHdr();
		ParameterNameValue paramType = hdr.findNameValParam(name);
		if (paramType != null) {
			result = paramType.getValueInt32();
		}
		return result;
	}

	protected void setIntToGenericHdr(String name, int value) {
		ParameterNameValue paramType = new ParameterNameValue();
		paramType.setName(name);
		paramType.setValueInt32(value);
		GenericDataHeader hdr = genericData.getHeader().getGenericDataHdr();
		hdr.addNameValParam(paramType);
	}

	protected void addParams(String type, List<ParameterNameValue> params) {
		GenericDataHeader hdr = genericData.getHeader().getGenericDataHdr();
		Iterator<ParameterNameValue> it = params.iterator();
		while (it.hasNext()) {
			ParameterNameValue param = it.next();
			param.setName(type + param.getName());
			hdr.addNameValParam(param);
		}
	}

}
