# makeFromNPUTE.pl, by Richie Lee (chichial@usc.edu), 1/25/2008
# ------------------------------------------------------------
# 
####################
# sub Functions
####################

sub setup{
  print "\nStarting the Perl program outputFromNPUTE.pl......\n\n";

  print STDOUT "Please enter a fileName file: ";
  chomp($FN = <STDIN>);
  print "Step A: Loading the source file '$FN'......\n\n";    
  open(SF,"$FN") or die "Can't open file\n";       

  print STDOUT "Please enter the main output file name: ";
  chomp($outputFile = <STDIN>);
  print "Step B: Locating the output file '$outputFile'.......\n\n";
  $start = (time)[0];   
}

sub closeFiles{
  $over = (time)[0];
  $totalSec = $over - $start;
  $totalMin = $totalSec / 60;
  print "\n***Total running time are $totalMin minutes***\n\n\n";  
  close(SF);close(SF1);close(SF2);close(SF3);close(SF4);close(SF5);    
}

sub makeFile{     
  $counter=0;
  $output1 = "output_ch1.txt";  
  $output2 = "output_ch2.txt";  
  $output3 = "output_ch3.txt";  
  $output4 = "output_ch4.txt";  
  $output5 = "output_ch5.txt";  
  open(SF1,"$output1") or die "Can't open output file\n";   
  open(SF2,"$output2") or die "Can't open output file\n";   
  open(SF3,"$output3") or die "Can't open output file\n"; 
  open(SF4,"$output4") or die "Can't open output file\n"; 
  open(SF5,"$output5") or die "Can't open output file\n";     
  @ch15=();
  @c15=();
  $c1=0;   
  while($line = <SF1>){
    $c1++;
    chomp($line);                  
    @c15 = split /,/, $line; 
    $s15 = scalar(@c15);
    for($i=0;$i<$s15;$i++){
    	$ch15[$i] .= $c15[$i].", ";    	
    }              	  	            
  }
  @c15=();  
  $c2=0;
  while($line = <SF2>){
    $c2++;
    chomp($line);                  
    @c15 = split /,/, $line; 
    $s15 = scalar(@c15);
    for($i=0;$i<$s15;$i++){
    	$ch15[$i] .= $c15[$i].", ";   
    }              	  	      	  	            
  } 
  @c15=();
  $c3=0;
  while($line = <SF3>){
    $c3++;
    chomp($line);                  
    @c15 = split /,/, $line; 
    $s15 = scalar(@c15);
    for($i=0;$i<$s15;$i++){
    	$ch15[$i] .= $c15[$i].", ";   
    }              	  	      	  	            
  }
  @c15=();
  $c4=0;
  while($line = <SF4>){
    $c4++;
    chomp($line);                  
    @c15 = split /,/, $line; 
    $s15 = scalar(@c15);
    for($i=0;$i<$s15;$i++){
    	$ch15[$i] .= $c15[$i].", ";   
    }              	  	      	  	            
  }
  @c15=();
  $c5=0;
  @ch5=();
  while($line = <SF5>){
    $c5++;
    chomp($line);                  
    @c15 = split /,/, $line; 
    $s15 = scalar(@c15);
    for($i=0;$i<$s15;$i++){    	  
    	$ch15[$i] .= $c15[$i].", ";   
    }              	  	      	  	            
  }  
  while($line = <SF>){  	            
      $counter++;
      chomp($line);                  
      @fileNames = split /\./, $line;     
      $fn1=$fileNames[0]."_base-calls.txt"; 
      $fn2=$fileNames[0]."_npute.txt";                     
      print "[$counter] Loading the file '$fn1'......\n";    
      open(SF0,"$fn1") or die "Can't open file\n";       
      open(result,">$fn2") or die "Can't open output file\n";   
      $counter1=0;
      @call = split /,/, $ch15[$counter-1];
      $s15 = scalar(@call);      
      $a=0;
      while($line1 = <SF0>){      	
      	if($counter1>0){
      	  chomp($line1);
      	  @data = split /\t/, $line1;
      	  print result "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t$data[7]\t$call[$a++]\n"; 
        }
        else{
          print result "$line1";
        }        
        $counter1++;  	                  	
      }      
      close(SF0);
      close(result);      
  }
}

sub printOriArray{
  $s = scalar(@_);          
  for($i=0;$i<$s;$i++){     
    print "\n$_[$i]";                           
  }
  print "\n";
}

####################
# Main Program
####################
setup;
makeFile;
closeFiles;

