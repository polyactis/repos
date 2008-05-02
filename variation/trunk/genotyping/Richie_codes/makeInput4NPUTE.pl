# makeNPUTE.pl, by Richie Lee (chichial@usc.edu), 1/25/2008
# ------------------------------------------------------------
# 
####################
# sub Functions
####################

sub setup{
  print "\nStarting the Perl program makeInput4NPUTE.pl......\n\n";

  print STDOUT "Please enter a source file name: ";
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
  close(SF);  
  close(result1);close(result2);close(result3);close(result4);close(result5);
}

sub makeFile{     
  $counter=0;
  $output1 = "ch1_".$outputFile;  
  $output2 = "ch2_".$outputFile;  
  $output3 = "ch3_".$outputFile;  
  $output4 = "ch4_".$outputFile;  
  $output5 = "ch5_".$outputFile;  
  open(result1,">$output1") or die "Can't open output file\n";   
  open(result2,">$output2") or die "Can't open output file\n";   
  open(result3,">$output3") or die "Can't open output file\n"; 
  open(result4,">$output4") or die "Can't open output file\n"; 
  open(result5,">$output5") or die "Can't open output file\n";     
  @ch1=();@ch2=();@ch3=();@ch4=();@ch5=();    
  while($line = <SF>){  	            
      $counter++;
      chomp($line);                  
      @fileNames = split /\./, $line;     
      $fn1=$fileNames[0]."_base-calls.txt";               
      print "[$counter] Loading the file '$fn1'......\n";    
      open(SF1,"$fn1") or die "Can't open file\n";       
      $counter1=0;
      $i1=0;$i2=0;$i3=0;$i4=0;$i5=0;
      while($line1 = <SF1>){
      	if($counter1>0){
      	  chomp($line1);
      	  @data = split /\t/, $line1;
      	  if($data[0] eq "1"){      	  
      	    $ch1[$i1] .= $data[8].", ";
      	    $i1++;
      	  }
      	  elsif($data[0] eq "2"){      	  
      	    $ch2[$i2] .= $data[8].", ";
      	    $i2++;
      	  }
      	  elsif($data[0] eq "3"){      	  
      	    $ch3[$i3] .= $data[8].", ";
      	    $i3++;
      	  }
      	  elsif($data[0] eq "4"){      	  
      	    $ch4[$i4] .= $data[8].", ";
      	    $i4++;
      	  }
      	  elsif($data[0] eq "5"){      	  
      	    $ch5[$i5] .= $data[8].", ";
      	    $i5++;
      	  } 
        }
      	$counter1++;  	                  	
      }
      close(SF1);      
  }
  $s1 = scalar(@ch1);for($i=0;$i<$s1;$i++){print result1 "$ch1[$i]\n";}
  $s2 = scalar(@ch2);for($i=0;$i<$s2;$i++){print result2 "$ch2[$i]\n";}
  $s3 = scalar(@ch3);for($i=0;$i<$s3;$i++){print result3 "$ch3[$i]\n";}
  $s4 = scalar(@ch4);for($i=0;$i<$s4;$i++){print result4 "$ch4[$i]\n";}
  $s5 = scalar(@ch5);for($i=0;$i<$s5;$i++){print result5 "$ch5[$i]\n";}
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

