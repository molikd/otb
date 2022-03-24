#!/bin/bash

error () { printf "%b\n" "[$(date)]: \e[91m$*\e[0m" >&2; exit 1; }
warn () { printf "%b\n" "[$(date)]: \e[93m$*\e[0m" >&2; }
state () { printf "%b\n" "[$(date)]: $*" 2>&1; }
describe () { printf "%b\n" "$*" 2>&1; }
pizzaz () { printf "%b\n" "\e[92m$*\e[0m" 2>&1; }
collect () { printf "%b\n" "\e[96m$*\e[0m" 2>&1; read var; return $var; }
version() { describe "otb: Only The Best (Genome Assemblies): v0.2.0"; exit 0;}

stop_check() {
  if [ -z "$force" ]; then
    printf "%b\n" "\e[96m$* (y/n)?:\e[0m" 2>&1;
    read cont
    case "$cont" in
      Y|y|yes|Yes|YES|YEs|yES|yEs|yeS )
      ;;
      N|n|no|No|NO|nO )
        error "\"NO\" entered, exiting";
      ;;
      sure )
        warn "sure, I mean, whatever..."
        stop_check
      ;;
      *)
        warn "Invalid input You entered: $cont";
	      stop_check
      ;;
    esac
  fi
}

display_header(){
printf "%b\n" "\e[7;40m ------------------------------------------------------------------------------------- \e[0m"
printf "%b\n" "\e[7;40m                                              ,,                                       \e[0m"
printf "%b\n" "\e[7;40m                                       mm    *MM                                       \e[0m"
printf "%b\n" "\e[7;40m                                       MM     MM                                       \e[0m"
printf "%b\n" "\e[7;40m                            ,pW\"Wq.  mmMMmm   MM,dMMb.                                 \e[0m"
printf "%b\n" "\e[7;40m                           6W'   \`Wb   MM     MM    \`Mb                                \e[0m"
printf "%b\n" "\e[7;40m                           8M     M8   MM     MM     M8                                \e[0m"
printf "%b\n" "\e[7;40m                           YA.   ,A9   MM     MM.   ,M9                                \e[0m"
printf "%b\n" "\e[7;40m                            \`Ybmd9'    \`Mbmo  P^YbmdP'                                 \e[0m"
printf "%b\n" "\e[7;40m                                                                                       \e[0m"
printf "%b\n" "\e[7;40m ------------------------------------------------------------------------------------- \e[0m"
printf "%b\n" "\e[7;40m               _   ._   |        _|_  |_    _     |_    _    _  _|_  |                 \e[0m"
printf "%b\n" "\e[7;40m              (_)  | |  |  \/     |_  | |  (/_    |_)  (/_  _>   |_  o                 \e[0m"
printf "%b\n" "\e[7;40m                           /                                                           \e[0m"
printf "%b\n" "\e[7;40m  /   _    _   ._    _   ._ _    _      _.   _   _   _   ._ _   |_   |  o   _    _  \  \e[0m"
printf "%b\n" "\e[7;40m |   (_|  (/_  | |  (_)  | | |  (/_    (_|  _>  _>  (/_  | | |  |_)  |  |  (/_  _>   | \e[0m"
printf "%b\n" "\e[7;40m  \   _|                                                                            /  \e[0m"
printf "%b\n" "\e[7;40m ------------------------------------------------------------------------------------- \e[0m"
}
