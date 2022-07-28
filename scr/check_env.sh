#!/bin/bash

if [ command -v io.sh &> /dev/null ] || [ -f "io.sh" ]; then
  source io.sh
elif [ -f "scr/io.sh" ]; then
  source scr/io.sh
else
  echo >&2 "[$(date)] I require io.sh to run, but io.sh has not been found, aborting"; exit 1;
fi

state "checking computing environment"

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
printf "%b\n" "\e[96m[SINGULARITY_LOCALCADHEDIR]\e[0m Singularity local cache directory:" 2>&1;
pizzaz "$SINGULARITY_LOCALCADHEDIR"
printf "%b\n" "\e[96m[SINGULARITY_CACHEDIR]\e[0m Singularity cache directory:" 2>&1;
pizzaz "$SINGULARITY_CACHEDIR"
printf "%b\n" "\e[96m[SINGULARITY_TMPDIR]\e[0m Singularity temporary directory:" 2>&1;
pizzaz "$SINGULARITY_TMPDIR"
printf "%b\n" "\e[96m[NXF_SINGULARITY_CACHEDIR]\e[0m Nextflow's Singularity cache directory:" 2>&1;
pizzaz "$NXF_SINGULARITY_CACHEDIR"
printf "%b\n" "\e[96m[NXF_SINGULARITY_LIBRARYDIR]\e[0m Nextflow's Singularity Library:" 2>&1;
pizzaz "$NXF_SINGULARITY_LIBRARYDIR"
describe ""
state "Singularity Things"
printf "%b\n" "\e[96m[which singularity]\e[0m Which Singularity?" 2>&1;
pizzaz "$( which singularity )"
printf "%b\n" "\e[96m[singularity version]\e[0m Singularity versioning:" 2>&1;
pizzaz "$( singularity version )"
describe ""
state "Nextflow Things"
printf "%b\n" "\e[96m[which nextflow]\e[0m Which Nextflow?" 2>&1;
pizzaz "$( which nextflow )"
printf "%b\n" "\e[96m[nextflow info]\e[0m Nextflow versioning:" 2>&1;
pizzaz "$( nextflow info )"
