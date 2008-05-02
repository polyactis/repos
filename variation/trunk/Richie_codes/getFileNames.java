import java.io.*;
public class getFileNames{
    static public void main (String[] args)
    {
    	try{
    	String Path = "D:\\Program Files\\R\\R-2.6.0\\";
    	String dirName = "JB-atSNPtile-24S-3-12-08";
    	File dir = new File(Path+dirName);
    	FileWriter fw = new FileWriter(Path+dirName+"_fileNames.txt");
    	String[] fileNames;
    	
    	if(dir.isDirectory()){    	  
    	  fileNames = dir.list();
    	  for(int i=0; i<fileNames.length; i++){
    	    fw.write(fileNames[i] + "\n");
    	  }
    	  fw.close();
    	}   
        }
        catch(Exception e){
          System.out.println("Exception: " + e);
        }     
    }
}
