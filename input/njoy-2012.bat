ECHO OFF
SET JOB=%1
SET ISOTOPE=%2
SET SYSTEMTYPE=dos
SET NJOY_DIR=\sync_Projects-Git\NJOY-2012\snl-work
cd \sync_Projects-Git\NJOY-2012\snl-work
perl \sync_Projects-Git\NJOY-2012\input\njoy-2012.prl %1 
rem if errorlevel 1 goto :Error
cd \sync_Projects-Git\NJOY-2012\input
ECHO Normal_NJOY-2012_batch_file_termination-Bye!
goto :END
rem :Error
rem echo NJOY-2012 batch file error detected
:END
