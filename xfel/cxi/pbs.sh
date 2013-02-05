#! /bin/sh

# $Id$

# Path to the chosen pyana script.  This should not need to be
# changed.  According to Marc Messerschmidt following ana-current
# should always be fine, unless one really wants to make sure
# everything is kept at the point where one started developing.  Do
# not use the shell's built-in which(1), which may give a relative
# path.
PYANA=`/usr/bin/which cxi.pyana`
if ! test -x "${PYANA}"; then
    echo "Cannot execute ${PYANA}" > /dev/stderr
    exit 1
fi

cleanup_and_exit() {
    exit ${1}
}
trap "cleanup_and_exit 1" HUP INT QUIT TERM

args=`getopt c:o:p:q:r:x: $*`
if test $? -ne 0; then
    echo "Usage: pbs.sh -c config -r run-num [-o output] [-p num-cpu] [-q queue] [-x exp]" > /dev/stderr
    cleanup_and_exit 1
fi

set -- ${args}
while test ${#} -ge 0; do
    case "${1}" in
        -c)
            cfg="${2}"
            if ! test -r "${cfg}" 2> /dev/null; then
                echo "config must be a readable file" > /dev/stderr
                cleanup_and_exit 1
            fi
            shift
            shift
            ;;

        -o)
            out=`readlink -fn "${2}"`
            if test -e "${out}" -a ! -d "${out}" 2> /dev/null; then
                echo "output exists but is not a directory" > /dev/stderr
                cleanup_and_exit 1
            fi
            test -d "${out}" 2> /dev/null || \
                echo "output directory will be created" > /dev/stderr
            shift
            shift
            ;;

        -p)
            if ! test "${2}" -gt 0 2> /dev/null; then
                echo "num-cpu must be positive integer" > /dev/stderr
                cleanup_and_exit 1
            fi
            nproc="${2}"
            shift
            shift
            ;;

        -q)
            queue="${2}"
            shift
            shift
            ;;

        -r)
            # Set ${run} to a zero-padded, four-digit string
            # representation of the integer.
            if ! test "${2}" -gt 0 2> /dev/null; then
                echo "run-num must be positive integer" > /dev/stderr
                cleanup_and_exit 1
            fi
            run=`echo "${2}" | awk '{ printf("%04d", $1); }'`
            shift
            shift
            ;;

        -x)
            exp="${2}"
            shift
            shift
            ;;

        --)
            shift
            break
            ;;
    esac
done

# Ensure the two mandatory arguments given, and no extraneous
# arguments are present.  XXX Since the corresponding options are not
# optional, they should perhaps be positional arguments instead?
if test -z "${cfg}" -o -z "${run}"; then
    echo "Must specify -c and -r options" > /dev/stderr
    cleanup_and_exit 1
fi
if test "${#}" -gt 0; then
    echo "Extraneous arguments" > /dev/stderr
    cleanup_and_exit 1
fi

# Take ${exp} from the environment unless overridden on the command
# line, and find its absolute path.
test -n "${EXP}" -a -z "${exp}" && exp="${EXP}"
exp=`find "/global/project/projectdirs/lcls/CXI" -maxdepth 2 -name "${exp}"`
if ! test -d "${exp}" 2> /dev/null; then
    echo "Could not find experiment subdirectory for ${exp}" > /dev/stderr
    cleanup_and_exit 1
fi

# Construct an absolute path to the directory with the XTC files as
# well as a sorted list of unique, comma-separated stream numbers for
# ${run}.  XXX May need some filtering as suggested by Amedeo Perazzo.
xtc="${exp}"
streams=`ls "${xtc}"/e*-r${run}-s* 2> /dev/null \
    | sed -e "s:.*-s\([[:digit:]]\+\)-c.*:\1:"  \
    | sort -u                                   \
    | awk '{ sep = NR > 1 ? "," : ""; printf("%s%d", sep, $1); }'`
if test -z "${streams}"; then
    echo "No streams in ${xtc}" > /dev/stderr
    cleanup_and_exit 1
fi

# If ${nproc} is not given on the the command line, fall back on
# num-cpu from ${cfg}.  Otherwise, the number of processes per host
# should be between 7 and 9 according to Marc Messerschmidt.  Using
# only two processors may decrease performance, because distributing
# data from the master process to a single worker process introduces
# overhead.
if test -z "${nproc}"; then
    nproc=`awk -F= '/^[[:space:]]*num-cpu[[:space:]]*=/ { \
                        printf("%d\n", $2);               \
                    }' "${cfg}"`
    test "${nproc}" -gt 0 2> /dev/null || nproc="7"
fi
if ! test ${nproc} != 2 2> /dev/null; then
    echo "Warning: running with two processors makes no sense" > /dev/stderr
fi

# If no queue is given on the command line then submit to default
# queue.
test -z "${queue}" && queue="regular"

# Unless specified on the command line, set up the output directory as
# a subdirectory named "results" within the experiment's scratch
# space.  All actual output will be written to the next available
# three-digit trial directory for the run.
if test -z "${out}"; then
    out="${exp}/scratch/results"
fi
out="${out}/r${run}"
trial=`mkdir -p "${out}" ;                 \
     find "${out}" -maxdepth 1             \
                   -name "[0-9][0-9][0-9]" \
                   -printf "%f\n" |        \
     sort -n | tail -n 1`
if test -z "${trial}"; then
    trial="000"
else
    if test "${trial}" -eq "999"; then
        echo "Error: Trial numbers exhausted" > /dev/stderr
        cleanup_and_exit 1
    fi
    trial=`expr "${trial}" \+ 1 | awk '{ printf("%03d", $1); }'`
fi
out="${out}/${trial}"

# Write a configuration file for the analysis of each stream by
# substituting the directory names with appropriate directories in
# ${out}, and appending the stream number to the base name.  XXX Dump
# the environment in here, too?  XXX What about an option to submit
# all streams to a single host, so as to do averaging?
mkdir -p "${out}"
for s in `echo "${streams}" \
    | tr -s ',' '\n'        \
    | awk -F, '{ printf("%02d ", $1); }'`; do
    sed -e "s:\([[:alnum:]]\+\)\(_dirname[[:space:]]*=\).*:\1\2 ${out}/\1:"    \
        -e "s:\([[:alnum:]]\+_basename[[:space:]]*=.*\)[[:space:]]*:\1s${s}-:" \
        "${cfg}" > "${out}/pyana_s${s}.cfg"
done

# Create all directories for the output from the analysis.  This
# eliminates a race condition when run in parallel.
mkdir -p "${out}/stdout"
awk -F=                                    \
    '/^[[:space:]]*[[:alnum:]]+_dirname[[:space:]]*=/ { \
         gsub(/^ /, "", $2);                            \
         gsub(/ $/, "", $2);                            \
         printf("\"%s\"\n", $2);                        \
     }' "${out}"/pyana_s[0-9][0-9].cfg | sort -u | xargs mkdir -p

# Because the streams are submitted for processing as a job array with
# one stream per host and ${nproc} processors per node, there are no
# individual submit files, like in lsf.sh.  The PBS script is not to
# be executed directly, but has to be passed to qsub.  XXX PBS seems
# to create output files stdout-[0-9]+; can that be changed to
# s00.out, s01.out, etc, for conformance with lsf.sh?  XXX The maximum
# wallclock limit is hardcoded to suit the carver regular queue.
#
# XXX From lsf.sh: allocate no more than ${nproc} processors.  Allow
# the job to start if at least one processor is available on the host.
cat > "${out}/submit.pbs" << EOF
#PBS -N r${run}
#PBS -j oe
#PBS -l gres=project:1,nodes=1:ppn=${nproc},walltime=04:00:00
#PBS -m a
#PBS -o ${out}/stdout/stdout
#PBS -q ${queue}

cd "\${PBS_O_WORKDIR}"

export SIT_ARCH="x86_64-gentoo-gcc463-opt"
export SIT_DATA="/global/project/projectdirs/lcls/psdm/data"
export SIT_RELEASE="ana-current"
export SIT_ROOT="/global/project/projectdirs/lcls/psdm-hattne"

export LD_LIBRARY_PATH=\${SIT_ROOT}/arch/x86_64-gentoo-gcc463-opt/lib:\${LD_LIBRARY_PATH}
export PYTHONPATH=\${SIT_ROOT}/arch/x86_64-gentoo-gcc463-opt/python:\${PYTHONPATH}

stream=\`echo "\${PBS_ARRAYID}" | awk '{ printf("%02d", \$1); }'\`

"${PYANA}" \
  -c "${out}/pyana_s\${stream}.cfg" \
  -p "\${PBS_NUM_PPN}" \
  "${xtc}"/e*-r${run}-s\${stream}-c*
EOF

# Copy the configuration files and the submission script to ${out}.
# Submit the analysis of all streams to the queueing system from
# ${NODE}.
qsub -t "${streams}" "${out}/submit.pbs"

echo "Output directory: ${out}"
cleanup_and_exit 0
