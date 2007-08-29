<?php

if ($_POST['_submit_check']) 
{
	$username="root";
	$password="";
	$database="retreat";
	$link=mysqli_connect(localhost,$username,$password,$database);
	if(!$link){ echo " mysql connection not working";}
	$name=$_POST['name'];
	echo  "<script  type='text/javascript'>
           alert('$name : Your abstract has been submited. Thank you');
          </script>"; 
	$email=$_POST['email'];
	$pi=$_POST['pi'];
	$Pref=$_POST['Pref'];
	$title=$_POST['title'];
	$author=$_POST['author_list'];
	$abs=$_POST['abs'];
	$query=" insert INTO abstract (name,email,pi,pref,title,author_list,abstract) VALUES('$name','$email','$pi','$Pref','$title','$author','$abs')";  
	$result=mysqli_query($link,$query);
	if(!result){ echo "error in query";}
	
	
	$strLocation = "./retreat.html";
	header("location:".$strLocation);
	mysql_close(); 
}
?>

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<title>MCB Retreat</title>
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
<form method="post" name="abstract" action="<?php echo $_SERVER['PHP_SELF'];?>">
<table width="800" height="747" border="0" align="center" cellpadding="0" cellspacing="0">
  <tr>
    <td height="747" align="center" valign="top" background="scroll_back.jpg"><table width="800" height="572" border="0" align="center" cellpadding="0" cellspacing="0">
      <tr>
        <td width="162" height="96">&nbsp;</td>
        <td width="116">&nbsp;</td>
        <td width="93">&nbsp;</td>
        <td width="343">&nbsp;</td>
        <td width="86">&nbsp;</td>
      </tr>
      <tr>
        <td height="39" colspan="5" align="center" valign="middle"><h1><span class="style3">MCB Retreat 2006</span> </h1></td>
        </tr>
      
      <tr>
        <td height="30" colspan="5" align="center" valign="middle"><p class="style8">Abstract</p></td>
        </tr>
      <tr>
        <td height="34">&nbsp;</td>
        <td colspan="2" align="left"><strong>Name:</strong></td>
        <td align="left"> <input type="text" name="name" /></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
        <td height="32">&nbsp;</td>
        <td colspan="2" align="left"><strong>E-mail:</strong></td>
        <td align="left"><input type="text" name="email" /></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
        <td height="32">&nbsp;</td>
        <td colspan="2" align="left"><strong>PI's Name: </strong></td>
        <td align="left"><input type="text" name="pi" /></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
        <td height="63" colspan="5"><div class="style7" style="margin:10px 40px 10px 60px">We have a limited number  of presentation times, so not all presentations can be accommodated.&nbsp;  Therefore some submitted presentations will have to be posters.&nbsp; We  will confirm with you by e-mail. </div></td>
        </tr>
     
	  <tr>
        <td height="33">&nbsp;</td>
        <td colspan="2" align="left"><strong>Preference: </strong></td>
        <td align="left"><select name="Pref"  >
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
        <td height="27" colspan="5" align="center" valign="middle"><input type="hidden" name="_submit_check" value="1"/><input type="submit" name="submit" ></td>
        </tr>
    </table></td>
  </tr>
</table>
</form>
</body>
</html>


