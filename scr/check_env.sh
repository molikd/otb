#!/bin/bash

if [ command -v io.sh &> /dev/null ] || [ -f "io.sh" ]; then
  source io.sh
elif [ -f "scr/io.sh" ]; then 
  source scr/io.sh
else
  echo >&2 "[$(date)] I require io.sh to run, but io.sh has not been found, aborting"; exit 1; 
fi

state "YOUR ENVIRONMENT"

describe "$( 
# thanks to linux mafia for this script http://linuxmafia.com/faq/Admin/release-files.html
OS=`uname -s`
REV=`uname -r`
MACH=`uname -m`

GetVersionFromFile()
{
    VERSION=`cat $1 | tr "\n" ' ' | sed s/.*VERSION.*=\ // `
}

if [ "${OS}" = "SunOS" ] ; then
    OS=Solaris
    ARCH=`uname -p` 
    OSSTR="${OS} ${REV}(${ARCH} `uname -v`)"
elif [ "${OS}" = "AIX" ] ; then
    OSSTR="${OS} `oslevel` (`oslevel -r`)"
elif [ "${OS}" = "Linux" ] ; then
    KERNEL=`uname -r`
    if [ -f /etc/redhat-release ] ; then
        DIST='RedHat'
        PSUEDONAME=`cat /etc/redhat-release | sed s/.*\(// | sed s/\)//`
        REV=`cat /etc/redhat-release | sed s/.*release\ // | sed s/\ .*//`
    elif [ -f /etc/SuSE-release ] ; then
        DIST=`cat /etc/SuSE-release | tr "\n" ' '| sed s/VERSION.*//`
        REV=`cat /etc/SuSE-release | tr "\n" ' ' | sed s/.*=\ //`
    elif [ -f /etc/mandrake-release ] ; then
        DIST='Mandrake'
        PSUEDONAME=`cat /etc/mandrake-release | sed s/.*\(// | sed s/\)//`
        REV=`cat /etc/mandrake-release | sed s/.*release\ // | sed s/\ .*//`
    elif [ -f /etc/debian_version ] ; then
        DIST="Debian `cat /etc/debian_version`"
        REV=""

    fi
    if [ -f /etc/UnitedLinux-release ] ; then
        DIST="${DIST}[`cat /etc/UnitedLinux-release | tr "\n" ' ' | sed s/VERSION.*//`]"
    fi

    OSSTR="${OS} ${DIST} ${REV}(${PSUEDONAME} ${KERNEL} ${MACH})"

fi

echo ${OSSTR}
)"

state "ENV Variables"
describe "[SINGULARITY_LOCALCADHEDIR] Singularity local cache directory: $SINGULARITY_LOCALCADHEDIR"
describe "[SINGULARITY_CACHEDIR] Singularity cache directory: $SINGULARITY_CACHEDIR"
describe "[SINGULARITY_TMPDIR] Singularity temporary directory: $SINGULARITY_TMPDIR"
describe "[NXF_SINGULARITY_CACHEDIR] Nextflow's Singularity cache directory: $NXF_SINGULARITY_CACHEDIR"
describe ""
state "Singularity Things"
describe "[which singularity] Which Singularity? $( which singularity )"
describe "[singularity version] Singularity versioning: $( singularity version )" 
describe ""
describe"Nextflow Things"
describe "[which nextflow] Which Nextflow? $( which nextflow )"
describe "[nextflow info] Nextflow versioning:"
describe "$( nextflow info )"
