import java.io.*;
import java.util.*;


import org.biojava.bio.program.sax.*;
import org.biojava.bio.program.ssbind.*;
import org.biojava.bio.search.*;
import org.biojava.bio.seq.db.*;
import org.xml.sax.*;
import org.biojava.bio.*;

public class BlastParser {
  /**
   * args[0] is assumed to be the name of a Blast output file
   */
  public static void main(String[] args) {
    try {
      //get the Blast input as a Stream
      InputStream is = new FileInputStream(args[0]);

      //make a BlastLikeSAXParser
      BlastLikeSAXParser parser = new BlastLikeSAXParser();

      //make the SAX event adapter that will pass events to a Handler.
      SeqSimilarityAdapter adapter = new SeqSimilarityAdapter();

      //set the parsers SAX event adapter
      parser.setContentHandler(adapter);

      //The list to hold the SeqSimilaritySearchResults
      List results = new ArrayList();

      //create the SearchContentHandler that will build SeqSimilaritySearchResults
      //in the results List
      SearchContentHandler builder = new BlastLikeSearchBuilder(results,
          new DummySequenceDB("queries"), new DummySequenceDBInstallation());

      //register builder with adapter
      adapter.setSearchContentHandler(builder);

      //parse the file, after this the result List will be populated with
      //SeqSimilaritySearchResults

      parser.parse(new InputSource(is));
      formatResults(results);
    }
    catch (SAXException ex) {
      //XML problem
      ex.printStackTrace();
    }catch (IOException ex) {
      //IO problem, possibly file not found
      ex.printStackTrace();
    }
  }
}
