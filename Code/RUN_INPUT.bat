@echo off
set "CURRENT_DIR=%~dp0"
matlab -nosplash -nodesktop -r "cd('%CURRENT_DIR%'); Input_GUI;" 