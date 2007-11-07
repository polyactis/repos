<?php

if ($_POST['_submit_check']) 
{
	$link = pg_connect("host=localhost dbname=yhdb user=yh password=123456")
		or die('Could not connect: ' . pg_last_error());
	$name=$_POST['name'];
	$email=$_POST['email'];
	$phone=$_POST['phone'];
	$address=$_POST['address'];
	$ride_type=$_POST['ride_type'];
	$no_of_people=$_POST['no_of_people'];
	$severity=$_POST['severity'];
	if ($address=='' || $email=='' || $name=='' ||$no_of_people=='')
	{
		echo  "<script  type='text/javascript'>
           alert('ERROR. Dear $name, Name or Email or Address or how many people is empty.');
          </script>";
	}
	else
	{
	$query=" insert INTO retreat.carpool (name, email, phone, address, ride_type, no_of_people, severity) VALUES('$name','$email','$phone', '$address', '$ride_type', '$no_of_people','$severity')";
	$result=pg_query($link,$query) or die('Query failed: ' . pg_last_error());
	pg_close($link);
	echo  "<script  type='text/javascript'>
           alert('$name : Thank you for your submission!');
          </script>";
	}
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
<form method="post" name="carpool" action="<?php echo $_SERVER['PHP_SELF'];?>">
<table width="800" height="747" border="0" align="center" cellpadding="0" cellspacing="0">
  <tr>
    <td height="747" align="center" valign="top" background="scroll_back.jpg"><table width="800" height="535" border="0" align="center" cellpadding="0" cellspacing="0">
      <tr>
        <td width="161" height="10">&nbsp;</td>
        <td width="117">&nbsp;</td>
        <td width="93">&nbsp;</td>
        <td width="343">&nbsp;</td>
        <td width="86">&nbsp;</td>
      </tr>
      <tr>
        <td height="30" colspan="5" align="center" valign="middle"><h1><span class="style3">MCB Retreat 2007</span> </h1></td>
        </tr>
      <tr>
        <td height="30" colspan="5" align="center" valign="middle"><p class="style8">Carpool</p></td>
        </tr>
      <tr>
        <td height="34">&nbsp;</td>
        <td colspan="2" align="left"><strong>*Name:</strong></td>
        <td align="left"> <input type="text" name="name" /></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
        <td height="32">&nbsp;</td>
        <td colspan="2" align="left"><strong>*E-mail:</strong></td>
        <td align="left"><input type="text" name="email" /></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
        <td height="32">&nbsp;</td>
        <td colspan="2" align="left"><strong>Phone: </strong></td>
        <td align="left"><input type="text" name="phone" /></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
        <td height="32">&nbsp;</td>
        <td colspan="2" align="left"><strong>*Address(for coordination): </strong></td>
        <td align="left"><input type="text" name="address" /></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
	<td>&nbsp;</td>
        <td colspan="2" align="left"><span class="style3">Need/Provide ride: </span></td>
	<td align="left"><select name="ride_type">
          <option>Need</option>
          <option>Provide</option>
        </select></td>
	<td>&nbsp;</td>
      </tr>
      <tr>
        <td height="32">&nbsp;</td>
        <td colspan="2" align="left"><strong>*How Many People you need/provide: </strong></td>
        <td align="left"><input type="text" name="no_of_people" /></td>
        <td>&nbsp;</td>
      </tr>
      <tr>
	<td>&nbsp;</td>
        <td colspan="2" align="left"><strong>Severity of issue (only for people who need ride): </strong></td>
	<td align="left"><select name="severity">
          <option>No Car(have to)</option>
          <option>has car but prefer carpool</option>
        </select></td>
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
