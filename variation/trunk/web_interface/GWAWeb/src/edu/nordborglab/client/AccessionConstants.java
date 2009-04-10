package edu.nordborglab.client;

/**
 * Interface to represent the constants contained in resource bundle:
 * 	'/usr/local/home_ubuntu/crocea/script/variation/web_interface/GWAWeb/src/edu/nordborglab/client/AccessionConstants.properties'.
 */
public interface AccessionConstants extends com.google.gwt.i18n.client.Constants {
  
  /**
   * Translated "/Accession/findAccessions".
   * 
   * @return translated "/Accession/findAccessions"
   */
  @DefaultStringValue("/Accession/findAccessions")
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
   * Translated "".
   * 
   * @return translated ""
   */
  @DefaultStringValue("")
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
   * Translated "".
   * 
   * @return translated ""
   */
  @DefaultStringValue("")
  @Key("cwAccessionByNameName")
  String cwAccessionByNameName();
}
