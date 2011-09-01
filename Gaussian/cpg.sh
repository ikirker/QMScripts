#!/usr/bin/env bash
# Copies a Gaussian .com file, changing the checkpoint file spec %chk= to the equivalent of the name of the new file.
##

if [[ "$1" == "" || "$2" == "" || "$3" != "" ]]; then
  echo "Incorrect number of arguments: need a source file, and a destination filename." >&2
  exit 2
fi

RelativeFileInName=$1
RelativeFileOutName=$2

if [[ "$RelativeFileInName" == *[/]* ]]; then
  PathIn=${RelativeFileInName%/*}
else
  PathIn=./
fi

if [[ "$RelativeFileOutName" == *[/]* ]]; then
  PathOut=${RelativeFileOutName%/*}
else
  PathOut=./
fi

FileIn=${1##*/}
FileOut=${2##*/}

NameIn=${FileIn%\.com}
NameOut=${FileOut%\.com}

ExistingCheckPointLine=`grep '%chk=' $FileIn`
if [[ "${ExistingCheckPointLine%\.chk}" != "%chk=$NameIn" ]]; then
  echo "Warning: existing .com file $FileIn does not have a conventional checkpoint file name." >&2
fi

echo "Changed %chk line to:" >&2
sed -ne "s/%chk=$NameIn/%chk=$NameOut/
w $RelativeFileOutName
T end 
p 
:end" $RelativeFileInName >&2 || echo "Could not complete - check other errors." >&2


