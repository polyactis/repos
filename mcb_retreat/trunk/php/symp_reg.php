<?php

if ($_POST['_submit_check']) 
{
	$link = pg_connect("host=localhost dbname=yhdb user=yh password=123456")
		or die('Could not connect: ' . pg_last_error());
	$name=trim($_POST['name']);
	$email=$_POST['email'];
	$pi=$_POST['pi'];
	$Pref=$_POST['Pref'];
	$title=$_POST['title'];
	$author=$_POST['author_list'];
	$abs=$_POST['abs'];
	//2008-04-14 file is in the _FILES array
	$abs_file=$_FILES['abs_file']["tmp_name"];
	
	if ($email=='' || $name=='')
	{
		echo  "<script  type='text/javascript'> 
           alert('Eror: Name and E-mail are required.');
           </script>";
	}
	elseif ($Pref!='None' && ( $title=='' || $pi=='' || $author=='' || ($abs=='' && $abs_file=='')))
	{
		echo  "<script  type='text/javascript'> 
           alert('Eror: since your preference is not None, please fill in other fields.');
           </script>";
	}
	else
	{
		if ($abs == '')
		{	
			/*
			echo "<pre>";
			echo "filename: " . $_FILES['abs_file']['name'] . '\n';
			echo $abs_file . "\n";
			echo $_FILES['abs_file']['name'] . '\n';
			echo $_FILES['abs_file']['size'] . '\n';
			echo $_FILES['abs_file']["type"] . '\n';
			*/
			
			$abs_file_content_type = $_FILES['abs_file']["type"];
			/*
			switch ($abs_file_content_type) {
			case "application/pdf":
				$abs_file_extension='pdf';
			    break;
			case 'application/msword':
				$abs_file_extension='doc';
				break;
			case 'text/plain':
				$abs_file_extension='txt';
				break;
			case 'text/rtf':
				$abs_file_extension='rtf';
				break;
			default:
				echo  "<script  type='text/javascript'> 
				alert('Error: Dear $name, your file has the unsupported extension: $abs_file_extension.');
				</script>";
				die('');
			}
			*/
			
			//2008-04-14 easier way to get file extension. but not used.
			$path_parts = pathinfo($_FILES['abs_file']['name']);
			$abs_file_extension = strtolower($path_parts['extension']);
			//echo "File Extension: " . $abs_file_extension . '\n';
			
			$abs_binary = file_get_contents($abs_file) or die('Reading file ' . $abs_file . ' failed.');
			$abs_file_data = base64_encode($abs_binary);	//base64_decode() is for the reverse action to get the file data out.
			$abs_filename = $name . "." . $abs_file_extension;
			$abs_file_size = $_FILES['abs_file']["size"];
			$query="insert INTO retreat.abstract (name,email,pi,pref,title,author_list, abstract_file_data, abstract_file_content_type, abstract_filename, abstract_file_size) 
				VALUES('$name','$email','$pi','$Pref','$title','$author','$abs_file_data', '$abs_file_content_type', '$abs_filename', '$abs_file_size')";
			$result=pg_query($link,$query) or die('DB insertion failed: ' . pg_last_error());
			
			//2007-04-14 upload the file
			$target_path = "uploads/";
			$target_path = $target_path . $abs_filename;
			move_uploaded_file($abs_file, $target_path);
		}
		else
		{
			$query=" insert INTO retreat.abstract (name,email,pi,pref,title,author_list,abstract) VALUES('$name','$email','$pi','$Pref','$title','$author','$abs')";
			$result=pg_query($link,$query) or die('DB insertion failed: ' . pg_last_error());
		}
		pg_close($link);
		echo  "<script  type='text/javascript'>
			alert('$name : Your abstract has been submited. Thank you!');
			</script>";
	}
}
?>

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<title>Biology Interdepartmental Graduate Symposium</title>
<style type="text/css">
<!--
.style3 {
	font-family: "Monotype Corsiva";
	color: #FF0000;
}
.style7 {
	color: #000000;
	font-family: "Times New Roman", Times, serif;
	font-weight: bold;
}
.style8 {
	font-size: 16px;
	font-family: "Lucida Sans Unicode";
	font-weight: bold;
}
-->
</style>
</head>

<body>
<form enctype="multipart/form-data" method="post" name="abstract" action="<?php echo $_SERVER['PHP_SELF'];?>">
<table width="800" height="747" border="0" align="center" cellpadding="0" cellspacing="0">
  <tr>
    <td height="747" align="center" valign="top"><table width="800" height="572" border="0" align="center" cellpadding="0" cellspacing="0">
      <tr>
        <td width="162" height="80">&nbsp;</td>
        <td width="116">&nbsp;</td>
        <td width="93">&nbsp;</td>
        <td width="343">&nbsp;</td>
        <td width="86">&nbsp;</td>
      </tr>
      <tr>
        <td height="39" colspan="5" align="center" valign="middle"><h1><span class="style3">Biology Interdepartmental Graduate Symposium</span> </h1></td>
        </tr>
      
      <tr>
        <td height="34">&nbsp;</td>
        <td colspan="2" align="left"><strong>Name<span class="style3">*</span>:</strong></td>
        <td align="left"> <input type="text" name="name" /></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
        <td height="32">&nbsp;</td>
        <td colspan="2" align="left"><strong>E-mail<span class="style3">*</span>:</strong></td>
        <td align="left"><input type="text" name="email" /></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
        <td height="63" colspan="5"><div class="style7" style="margin:10px 40px 10px 60px">All presentations must be approved by your advisor. We have a limited number of presentation slots but ample space for posters. We especially encourage graduate students who have graduated recently or are planning to graduate soon to give a talk. If you plan to only attend, please select None below. We will confirm with you by e-mail. </div></td>
        </tr>
     
	  <tr>
        <td height="33">&nbsp;</td>
        <td colspan="2" align="left"><strong>Preference: </strong></td>
        <td align="left"><select name="Pref"  >
          <option>None</option>
          <option>Poster</option>
          <option>Presentation</option>
        </select></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
        <td height="33">&nbsp;</td>
        <td colspan="2" align="left"><strong>Title: </strong></td>
        <td align="left">
          
            <input name="title" type="text" size="55" />        </td>
        <td>&nbsp;</td>
      </tr>
      <tr>
        <td height="32">&nbsp;</td>
        <td colspan="2" align="left"><strong>PI's Name: </strong></td>
        <td align="left"><input type="text" name="pi" /></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
        <td height="33">&nbsp;</td>
        <td colspan="2" align="left"><strong>Author List:</strong></td>
        <td align="left"> <input name="author_list" type="text" size="55" /></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
        <td height="120">&nbsp;</td>
        <td colspan="2" align="left" valign="middle"><strong>Abstract: </strong></td>
        <td align="left"> <textarea name="abs" cols="52" rows="7"></textarea></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
        <td height="33">&nbsp;</td>
        <td colspan="2" align="left"><strong><span class="style3">OR</span> Abstract File: </strong></td>
        <td align="left"> <input TYPE="FILE" NAME="abs_file" SIZE="55"></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
        <td height="27" colspan="5" align="center" valign="middle"><input type="hidden" name="_submit_check" value="1"/><input type="submit" name="submit" ></td>
        </tr>
    </table></td>
  </tr>
</table>
</form>
</body>
</html>


