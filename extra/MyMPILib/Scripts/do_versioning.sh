#!/bin/bash
GITHASH=`git rev-parse HEAD`
GITSHORTHASH=`git rev-parse --short HEAD`
GITCHANGEDFILES=`git diff-index --name-only HEAD`

export GITSHORTHASH=$GITSHORTHASH

echo "character(len=*), parameter :: MyMPILib_Version = '${GITHASH}'" > ./Internal/version.f90
echo "Versioning MyMPILib..."

if [ -n "$GITCHANGEDFILES" ]; then
    echo 'character(len=*), parameter :: MyMPILib_Version_Additional = "WARNING, &' >> ./Internal/version.f90
    echo "&THERE ARE UNCOMMITTED CHANGES. Run may not be reproduceable: &" >> ./Internal/version.f90

    while read -r line; do
        echo "&${line} &" >> ./Internal/version.f90
    done <<< "$GITCHANGEDFILES"
    echo '&"' >> ./Internal/version.f90

else

    echo 'character(len=*), parameter :: MyMPILib_Version_Additional = ""' >> ./Internal/version.f90

fi
