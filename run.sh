#!/bin/bash

echo "Starting job"

if [[ "$HOSTNAME" == *"triolith"* ]]
then
  echo "We're on triolith!";

  p='1.85'
  mu='-0.011'
  a='1'
  N='100000000000'
  t='50'

  sbatch run_triolith.sh -p ${p} -mu ${mu} -a ${a} -N ${N} -t ${t}
else
  p='0'
  mu='-1.1'
  a='1'
  N='5000000000'
  t='30'

  while getopts 'p:m:a:N:t:' flag; do
    case "${flag}" in
      p) p="${OPTARG}" ;;
      m) mu="${OPTARG}" ;;
      a) a="${OPTARG}" ;;
      N) N="${OPTARG}" ;;
      t) t="${OPTARG}" ;;
      *) error "Unexpected option ${flag}" ;;
    esac
  done

  echo ./bin/run -p ${p} -mu ${mu} -a ${a} -N ${N} -t ${t}
  ./bin/run -p ${p} -mu ${mu} -a ${a} -N ${N} -t ${t}
fi