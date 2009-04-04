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
