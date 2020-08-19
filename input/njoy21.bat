ECHO OFF
SET JOB=%1
SET ISOTOPE=%2
SET SYSTEMTYPE=dos
SET NJOY_DIR=\sync_LINUX\NJOY21\snl-work
cd \sync_Linux\NJOY21\snl-work
perl \sync_Linux\NJOY21\input\njoy21.prl %1 
rem ./njoy21 <input >output
rem if errorlevel 1 goto :Error
cd \sync_Linux\NJOY21\input
ECHO Normal_NJOY21_batch_file_termination-Bye!
goto :END
rem :Error
rem echo NJOY21 batch file error detected
:END
