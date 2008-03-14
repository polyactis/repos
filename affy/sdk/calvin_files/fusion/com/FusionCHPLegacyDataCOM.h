////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////


#pragma once
#include "resource.h"       // main symbols

#include "affx_fusion_com.h"
#include "FusionCHPLegacyData.h"

// CFusionCHPLegacyDataCOM

class ATL_NO_VTABLE CFusionCHPLegacyDataCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFusionCHPLegacyDataCOM, &CLSID_FusionCHPLegacyData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFusionCHPLegacyData, &IID_IFusionCHPLegacyData, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFusionCHPLegacyDataCOM()
	{
		chp = NULL;
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FUSIONCHPLEGACYDATA)


BEGIN_COM_MAP(CFusionCHPLegacyDataCOM)
	COM_INTERFACE_ENTRY(IFusionCHPLegacyData)
	COM_INTERFACE_ENTRY(IDispatch)
	COM_INTERFACE_ENTRY(ISupportErrorInfo)
END_COM_MAP()

// ISupportsErrorInfo
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);


	DECLARE_PROTECT_FINAL_CONSTRUCT()

	HRESULT FinalConstruct()
	{
		return S_OK;
	}

	void FinalRelease()
	{
		delete chp;
        chp = NULL;
	}

private:
	affymetrix_fusion_io::FusionCHPLegacyData *chp;

public:
	STDMETHOD(get_FileTypeIdentifier)(BSTR* pVal);
public:
	STDMETHOD(get_FileTypeIdentifiers)(VARIANT* pVal);
public:
	STDMETHOD(FromBase)(IFusionCHPData * baseChp, VARIANT_BOOL* pVal);
public:
	STDMETHOD(GetExpressionResults)(int index, IFusionExpressionProbeSetResults * pVal);
public:
	STDMETHOD(GetGenotypingResults)(int index, IFusionGenotypeProbeSetResults * pVal);
public:
	STDMETHOD(GetUniversalResults)(int index, IFusionUniversalProbeSetResults * pVal);
public:
	STDMETHOD(GetReseqResults)(IFusionResequencingResults * pVal);
public:
	STDMETHOD(Clear)(void);
public:
	STDMETHOD(GetHeader)(IFusionCHPHeader ** pVal);
public:
	STDMETHOD(get_FileId)(BSTR* pVal);
public:
    STDMETHOD(GetProbeSetName)(int index, BSTR* pVal);
};

OBJECT_ENTRY_AUTO(__uuidof(FusionCHPLegacyData), CFusionCHPLegacyDataCOM)
