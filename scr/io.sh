#!/bin/bash

error () { printf "%b\n" "[$(date)]: \e[91m$*\e[0m" >&2; exit 1; }
warn () { printf "%b\n" "[$(date)]: \e[93m$*\e[0m" >&2; }
state () { printf "%b\n" "[$(date)]: $*" 2>&1; }
describe () { printf "%b\n" "$*" 2>&1; } 
pizzaz () { printf "%b\n" "\e[92m$*\e[0m" 2>&1; }
collect () { printf "%b\n" "\e[96m$*\e[0m" 2>&1; read var; return $var; }
version() { describe "otb: Only The Best (Genome Assemblies): v0.1.0"; exit 0;}

stop_check() {
  if [ -z "$force" ]; then
    printf "%b\n" "\e[96m$* (y/n)?:\e[0m" 2>&1;
    read cont
    case "$cont" in
      Y|y|yes|Yes|YES )
      ;;
      N|n|no|No|NO )
        error "\"NO\" entered, exiting";
      ;;
      *) 
        warn "Invalid input You entered: $cont";
	    stop_check
      ;;
    esac
  fi
}

display_header(){
  printf "%b\n" "\e[7;40m                                                                         \e[0m"
  printf "%b\n" "\e[7;40m       .----------------. .----------------. .----------------.          \e[0m"
  printf "%b\n" "\e[7;40m      | .--------------. | .--------------. | .--------------. |         \e[0m"
  printf "%b\n" "\e[7;40m      | |     ____     | | |  _________   | | |   ______     | |         \e[0m"
  printf "%b\n" "\e[7;40m      | |   .'    \`.   | | | |  _   _  |  | | |  |_   _ \    | |         \e[0m"
  printf "%b\n" "\e[7;40m      | |  /  .--.  \  | | | |_/ | | \_|  | | |    | |_) |   | |         \e[0m"
  printf "%b\n" "\e[7;40m      | |  | |    | |  | | |     | |      | | |    |  __'.   | |         \e[0m"
  printf "%b\n" "\e[7;40m      | |  \  \`--'  /  | | |    _| |_     | | |   _| |__) |  | |         \e[0m"
  printf "%b\n" "\e[7;40m      | |   \`.____.'   | | |   |_____|    | | |  |_______/   | |         \e[0m"
  printf "%b\n" "\e[7;40m      | |              | | |              | | |              | |         \e[0m"
  printf "%b\n" "\e[7;40m      | '--------------' | '--------------' | '--------------' |         \e[0m"
  printf "%b\n" "\e[7;40m       '----------------' '----------------' '----------------'          \e[0m"
  printf "%b\n" "\e[7;40m      ____       __       ________         ___         __    __          \e[0m"
  printf "%b\n" "\e[7;40m     / __ \___  / /_ __  /_  __/ /  ___   / _ )___ ___/ /_  / /          \e[0m"
  printf "%b\n" "\e[7;40m    / /_/ / _ \/ / // /   / / / _ \/ -_) / _  / -_|_-< __/ /_/           \e[0m"
  printf "%b\n" "\e[7;40m   _\____/_//_/_/\_, /   /_/ /_//_/\__/ /____/\__/___|__/_(_) ___        \e[0m"
  printf "%b\n" "\e[7;40m  / ___/__ ___  /___/__ _  ___   / _ | ___ ______ __ _  / /  / (_)__ ___ \e[0m"
  printf "%b\n" "\e[7;40m / (_ / -_) _ \/ _ \/  ' \/ -_) / __ |(_-<(_-< -_)  ' \/ _ \/ / / -_|_-< \e[0m"
  printf "%b\n" "\e[7;40m \___/\__/_//_/\___/_/_/_/\__/ /_/ |_/___/___|__/_/_/_/_.__/_/_/\__/___/ \e[0m"
  printf "%b\n" "\e[7;40m                                                                         \e[0m"
}


