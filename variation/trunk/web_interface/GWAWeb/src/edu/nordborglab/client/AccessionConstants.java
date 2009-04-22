package edu.nordborglab.client;

/**
 * Interface to represent the constants contained in resource bundle:
 * 	'/usr/local/home_ubuntu/crocea/script/variation/web_interface/GWAWeb/src/edu/nordborglab/client/AccessionConstants.properties'.
 */
public interface AccessionConstants extends com.google.gwt.i18n.client.Constants {
  
  /**
   * Translated "/Accession/findAccessionsByID?id=".
   * 
   * @return translated "/Accession/findAccessionsByID?id="
   */
  @DefaultStringValue("/Accession/findAccessionsByID?id=")
  @Key("AccessionByIDURL")
  String AccessionByIDURL();

  /**
   * Translated "/Accession/findAccessionsByName?name=".
   * 
   * @return translated "/Accession/findAccessionsByName?name="
   */
  @DefaultStringValue("/Accession/findAccessionsByName?name=")
  @Key("AccessionByNameURL")
  String AccessionByNameURL();

  /**
   * Translated "/Accession/autoComplete".
   * 
   * @return translated "/Accession/autoComplete"
   */
  @DefaultStringValue("/Accession/autoComplete")
  @Key("AccessionSuggestOracleURL")
  String AccessionSuggestOracleURL();

  /**
   * Translated "AccessionConstants".
   * 
   * @return translated "AccessionConstants"
   */
  @DefaultStringValue("AccessionConstants")
  @Key("ClassName")
  String ClassName();

  /**
   * Translated "/Phenotype/getPhenotypeIcon".
   * 
   * @return translated "/Phenotype/getPhenotypeIcon"
   */
  @DefaultStringValue("/Phenotype/getPhenotypeIcon")
  @Key("GetPhenotypeIconURL")
  String GetPhenotypeIconURL();

  /**
   * Translated "/Phenotype/getPhenotypeValue".
   * 
   * @return translated "/Phenotype/getPhenotypeValue"
   */
  @DefaultStringValue("/Phenotype/getPhenotypeValue")
  @Key("GetPhenotypeValueURL")
  String GetPhenotypeValueURL();

  /**
   * Translated "diaplay all".
   * 
   * @return translated "diaplay all"
   */
  @DefaultStringValue("diaplay all")
  @Key("MapWithPhenotypeDisplayOption1")
  String MapWithPhenotypeDisplayOption1();

  /**
   * Translated "diaplay only accessions with phenotype".
   * 
   * @return translated "diaplay only accessions with phenotype"
   */
  @DefaultStringValue("diaplay only accessions with phenotype")
  @Key("MapWithPhenotypeDisplayOption2")
  String MapWithPhenotypeDisplayOption2();

  /**
   * Translated "diaplay only accessions without phenotype".
   * 
   * @return translated "diaplay only accessions without phenotype"
   */
  @DefaultStringValue("diaplay only accessions without phenotype")
  @Key("MapWithPhenotypeDisplayOption3")
  String MapWithPhenotypeDisplayOption3();

  /**
   * Translated "diaplay all, diaplay only accessions with phenotype, diaplay only accessions without phenotype".
   * 
   * @return translated "diaplay all, diaplay only accessions with phenotype, diaplay only accessions without phenotype"
   */
  @DefaultStringValue("diaplay all, diaplay only accessions with phenotype, diaplay only accessions without phenotype")
  @Key("MapWithPhenotypeDisplayOptions")
  String MapWithPhenotypeDisplayOptions();

  /**
   * Translated "/Phenotype/getPhenotypeMethodLs".
   * 
   * @return translated "/Phenotype/getPhenotypeMethodLs"
   */
  @DefaultStringValue("/Phenotype/getPhenotypeMethodLs")
  @Key("MapWithPhenotypeGetPhenotypeMethodLsURL")
  String MapWithPhenotypeGetPhenotypeMethodLsURL();

  /**
   * Translated "call method id:".
   * 
   * @return translated "call method id:"
   */
  @DefaultStringValue("call method id:")
  @Key("callMethodLabel")
  String callMethodLabel();

  /**
   * Translated "".
   * 
   * @return translated ""
   */
  @DefaultStringValue("")
  @Key("cwAccession250kDescription")
  String cwAccession250kDescription();

  /**
   * Translated "".
   * 
   * @return translated ""
   */
  @DefaultStringValue("")
  @Key("cwAccession250kLabel")
  String cwAccession250kLabel();

  /**
   * Translated "250k Accessions ".
   * 
   * @return translated "250k Accessions "
   */
  @DefaultStringValue("250k Accessions ")
  @Key("cwAccession250kName")
  String cwAccession250kName();

  /**
   * Translated "Get ID from \"By Name\" Panel. ID has to be in number.".
   * 
   * @return translated "Get ID from \"By Name\" Panel. ID has to be in number."
   */
  @DefaultStringValue("Get ID from \"By Name\" Panel. ID has to be in number.")
  @Key("cwAccessionByIDDescription")
  String cwAccessionByIDDescription();

  /**
   * Translated "<b>Enter an ecotype ID: </b> ".
   * 
   * @return translated "<b>Enter an ecotype ID: </b> "
   */
  @DefaultStringValue("<b>Enter an ecotype ID: </b> ")
  @Key("cwAccessionByIDLabel")
  String cwAccessionByIDLabel();

  /**
   * Translated "By ID".
   * 
   * @return translated "By ID"
   */
  @DefaultStringValue("By ID")
  @Key("cwAccessionByIDName")
  String cwAccessionByIDName();

  /**
   * Translated "Support wild character like .* or ? etc. \"V.*r\"".
   * 
   * @return translated "Support wild character like .* or ? etc. \"V.*r\""
   */
  @DefaultStringValue("Support wild character like .* or ? etc. \"V.*r\"")
  @Key("cwAccessionByNameDescription")
  String cwAccessionByNameDescription();

  /**
   * Translated "<b>Enter an ecotype name: </b> ".
   * 
   * @return translated "<b>Enter an ecotype name: </b> "
   */
  @DefaultStringValue("<b>Enter an ecotype name: </b> ")
  @Key("cwAccessionByNameLabel")
  String cwAccessionByNameLabel();

  /**
   * Translated "By Name".
   * 
   * @return translated "By Name"
   */
  @DefaultStringValue("By Name")
  @Key("cwAccessionByNameName")
  String cwAccessionByNameName();

  /**
   * Translated "gene list type id: ".
   * 
   * @return translated "gene list type id: "
   */
  @DefaultStringValue("gene list type id: ")
  @Key("geneListTypeLabel")
  String geneListTypeLabel();

  /**
   * Translated "phenotype method id: ".
   * 
   * @return translated "phenotype method id: "
   */
  @DefaultStringValue("phenotype method id: ")
  @Key("phenotypeMethodLabel")
  String phenotypeMethodLabel();

  /**
   * Translated "Sort the SNP matrix by which Principal Component:  ".
   * 
   * @return translated "Sort the SNP matrix by which Principal Component:  "
   */
  @DefaultStringValue("Sort the SNP matrix by which Principal Component:  ")
  @Key("whichPCLabel")
  String whichPCLabel();
}
