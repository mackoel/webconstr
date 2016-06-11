<!DOCTYPE html>
<!--
To change this license header, choose License Headers in Project Properties.
To change this template file, choose Tools | Templates
and open the template in the editor.
-->
<html>
    <head>
        <meta charset="UTF-8">
        <style>
            .error {color: #FF0000;}
        </style>
        <title>
            DEEP Web Interface
        </title>
<?php

// define variables and set to empty values

$nameErr = $startErr = $endErr = $coderegErr = "";
$name = $start = $end = $codereg = "";
$codereg = "0";
if ($_SERVER["REQUEST_METHOD"] == "POST") {
   if (empty($_POST["name"])) {
     $nameErr = "Name is required";
   } else {
     $name = test_input($_POST["name"]);
     // check if name only contains letters and whitespace
     if (!preg_match("/^[a-zA-Z ]*$/",$name)) {
       $nameErr = "Only letters are allowed";
     }
   }
  
   if (empty($_POST["start"])) {
     $startErr = "Start is required";
   } else {
     $start = test_input($_POST["start"]);
     // check if e-mail address is well-formed
//     if (!filter_var($start, FILTER_VALIDATE_EMAIL)) {
//       $startErr = "Invalid email format";
//    }
   }
    
   if (empty($_POST["end"])) {
     $end = "";
     $endErr = "End is required";
   } else {
     $end = test_input($_POST["end"]);
     // check if URL address syntax is valid (this regular expression also allows dashes in the URL)
//     if (!preg_match("/\b(?:(?:https?|ftp):\/\/|www\.)[-a-z0-9+&@#\/%?=~_|!:,.;]*[-a-z0-9+&@#\/%=~_|]/i",$website)) {
//       $endErr = "end URL";
//     }
   }

   if (empty($_POST["codereg"])) {
     $coderegErr = "codereg is required";
   } else {
     $codereg = test_input($_POST["codereg"]);
   }
}

function test_input($data) {
   $data = trim($data);
   $data = stripslashes($data);
   $data = htmlspecialchars($data);
   return $data;
}
?>

<h2>Datafile</h2>
<p><span class="error">* required field.</span></p>
<form method="post" action="<?php echo htmlspecialchars($_SERVER["PHP_SELF"]);?>">
   <h3>Minimizer settings</h3>
   Name: <input type="text" name="name" value="<?php echo $name;?>">
   <span class="error">* <?php echo $nameErr;?></span>
   <br><br>
   Start: <input type="text" name="start" value="<?php echo $start;?>">
   <span class="error">* <?php echo $startErr;?></span>
   <br><br>
   End: <input type="text" name="end" value="<?php echo $end;?>">
   <span class="error"><?php echo $endErr;?></span>
   <br><br>
   Codereg:
   <input type="radio" name="codereg" <?php if (isset($codereg) && $codereg=="1") echo "checked";?>  value="1">1
   <input type="radio" name="codereg" <?php if (isset($codereg) && $codereg=="0") echo "checked";?>  value="0">0
   <span class="error">* <?php echo $coderegErr;?></span>
   <br><br>
   <input type="submit" name="submit" value="Submit">
</form>

<?php
echo "<h2>Your Input:</h2>";
echo $name;
echo "<br>";
echo $start;
echo "<br>";
echo $end;
echo "<br>";
echo $codereg;

echo '<pre>';

// Выводит весь результат шелл-команды "ls", и возвращает 
// последнюю строку вывода в переменной $last_line. Сохраняет код возврата
// шелл-команды в $retval.
// Если вы собираетесь передавать функции пользовательские данные, 
// используйте функции escapeshellarg() или escapeshellcmd() для того, 
// чтобы пользователи не смогли обмануть систему, запустив произвольную команду.
// Мы допускаем здесь произвольное количество аргументов умышленно.
//$command = './configure '.$_POST['configure_options'];

//$escaped_command = escapeshellcmd($command);

//$cmd = 'ls -al';

//$array_output = array(1, 1, 1, 1,  1, 8 => 1,  4 => 1, 19, 3 => 13);
//ob_start();
//$last_line = system($cmd, $retval);
//$last_line = exec($cmd, $array_output, $retval);
//$dat = ob_get_clean();
//$last_line = system('ls '.escapeshellarg($dir), $retval);
// Выводим дополнительную информацию
	$savep = getenv("PATH");        // save old value
	$newp = "/usr/local/bin:/home/other1/bin:/home/ann/bin";  // extra paths to add
	if ($savep) { $newp .= ":$savep"; }           // append old paths if any
	putenv("PATH=$newp");        // set new value
	$saveld = getenv("LD_LIBRARY_PATH");        // save old value
	$newld = "/usr/local/lib:/home/other1/lib:/home/ann/lib";  // extra paths to add
	if ($saveld) { $newld .= ":$saveld"; }           // append old paths if any
	putenv("LD_LIBRARY_PATH=$newld");        // set new value
//	$last_line = exec("./webconstr.R --check_constr=1 --check_dnase=1 --constr_target=gt --constr_local_start=1462 --constr_local_end=2674 -o dir2 2>/dev/null", $ret);        // do system command;
// mycommand is loaded using
// libs in the new path list
	$cmd = "./webconstr.R --check_constr=1 --check_dnase=1 --constr_target=".$name." --constr_local_start=".$start." --constr_local_end=".$end." --check_codereg=".$codereg." -o dir2 2>/dev/null";
//	$cmd = "echo 'dd' > dir2/n";
	$last_line = exec($cmd, $retval);
//	$last_line = system($cmd, $retval);
	putenv("PATH=$savep"); 
	putenv("LD_LIBRARY_PATH=$saveld");
	echo '<img src="' . $last_line . '" />';
echo '
</pre>
<hr />Команда: ' . $cmd . '
<hr />Последняя строка вывода: ' . $last_line . '
<hr />Код возврата: ' . $retval;
?>

</body>
</html>

