<html>
 <head>
  <title>PHP Page</title>
 </head>
 <body>
    <?php echo '<p>Image</p>'; 
	$savep = getenv("PATH");        // save old value
	$newp = "/usr/local/bin";  // extra paths to add
	if ($savep) { $newp .= ":$savep"; }           // append old paths if any
	putenv("PATH=$newp");        // set new value
	$saveld = getenv("LD_LIBRARY_PATH");        // save old value
	$newld = "/usr/local/lib";  // extra paths to add
	if ($saveld) { $newld .= ":$saveld"; }           // append old paths if any
	putenv("LD_LIBRARY_PATH=$newld");        // set new value
	$last_line = exec("./webconstr.R --check_constr=1 --check_dnase=1 --constr_target=gt --constr_local_start=1462 --constr_local_end=2674 -o dir2 2>/dev/null", $ret);        // do system command;
// mycommand is loaded using
// libs in the new path list
	putenv("PATH=$savep"); 
	putenv("LD_LIBRARY_PATH=$saveld");
	echo '<img src="' . $last_line . '" />';
?>
 </body>
</html>
