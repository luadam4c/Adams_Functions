#!/bin/bash
#
# Requires:
#   adams_commands.sh (for checkdir)
#   compile_script.m
#
# Used by:
#
# History: 
#   2018-01-31 - Modified from /usr/bin/template.sh
#   2018-02-05 - Changed folder name to minEASE_Linux_$version
#   2018-02-22 - Added local.settings to be copied
#   2018-02-28 - Changed folder path
#   2018-03-15 - Always compile with /usr/local/MATLAB/R2015b/bin/mcc for now
#   2020-08-20 - Now uses the newest version of MATLAB available

################################################################################

## Source functions from files
source adams_commands.sh

################################################################################

## Precautions
set -u                      # stop the script if a variable is not assigned
set -e                      # stop the script if an error is encountered

## Deal with arguments
# Check number of arguments
[[ $# -lt 1 ]] && echo USAGE: bash $0 version && exit 1

# Read in arguments
version=$1

################################################################################

# Create command for matlab
command=$(echo "addpath('/home/Matlab/Adams_Functions');" \
            "addpath('/home/Matlab/Downloaded_Functions');" \
            "compile_script('minEASE', 'ExtraFileNames', 'local.settings');" \
            "exit;")

# Run MATLAB to compile code
matlab -nodisplay -nosplash -r "${command}"

# Move compiled code to a folder named by the version
folder=/media/shareX/minEASE/minEASE_Linux_$version/
checkdir "$folder"
mv minEASE run_minEASE.sh readme.txt requiredMCRProducts.txt \
    mccExcludedFiles.log "$folder"
cp -p local.settings "$folder"

# Change permissions to be group writeable and group executable
chmodRgw "$folder"
chmodRgX "$folder"

## Return exit status upon success
exit 0

################################################################################
## OLD CODE:

# folder=/media/shareX/share/minEASE/minEASE_Linux_$version/

################################################################################