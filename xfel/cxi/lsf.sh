#! /bin/sh

# This script executes several commands over ssh.  It is probably a
# good idea to have an ssh-agent(1) running.  XXX Check again with
# notes to see that all this is sane.
#
# Note: A valid AFS token can be obtained by "kinit" followed by
# "aklog".  This avoids the "job being submitted without an AFS token"
# warning.
#
# $Id$

# This script must be run from the SIT directory, which contains the
# .sit_release file, so that the relative PYTHONPATH set by sit_setup
# is valid.  XXX Wouldn't it make sense to have
# /reg/g/psdm/etc/ana_env.sh set an absolute path?  Could find the
# user's release directory from .sit_release file and cd to it in the
# submit.sh script.  No, that's much too slow!
if ! relinfo > /dev/null 2>&1; then
    echo "Must run this script from the SIT release directory" > /dev/stderr
    exit 1
fi

# Absolute path to the cxi directory.
CXI="/reg/data/ana11/cxi"

# A random host that has the scratch directory mounted.  psexport is
# preferred over psanafeh, since the latter is not accessible from
# everywhere.
NODE="psexport.slac.stanford.edu"

# Path to the chosen pyana script.  This should not need to be
# changed.  According to Marc Messerschmidt following ana-current
# should always be fine, unless one really wants to make sure
# everything is kept at the point where one started developing.
PYANA="${SIT_RELDIR}/${SIT_RELEASE}/arch/${SIT_ARCH}/bin/pyana"
if ! test -x "${PYANA}"; then
    echo "Cannot execute ${PYANA}" > /dev/stderr
    exit 1
fi

# Absolute path to the python interpreter.  This should not need to be
# changed.  Do not use the shell's built-in which(1), which may give a
# relative path.
PYTHON=`/usr/bin/which libtbx.python`

args=`getopt c:o:p:r:x: $*`
if test $? -ne 0; then
    echo "Usage: lsf.sh -c config -r runno [-o output] [-p num-cpu] [-x exp]" > /dev/stderr
    exit 1
fi

set -- ${args}
while test $# -ge 0; do
    case "$1" in
        -c)
            if ! test -r "$2" 2> /dev/null; then
                echo "config must be a readable file" > /dev/stderr
                exit 1
            fi
            cfg="$2"
            shift
            shift
            ;;

        -o)
            if ssh ${NODE} "test -e \"$2\" -a ! -d \"$2\" 2> /dev/null"; then
                echo "output exists but is not a directory" > /dev/stderr
                exit 1
            fi
            ssh ${NODE} "test -d \"$2\" 2> /dev/null" || \
                echo "output directory will be created" > /dev/stderr
            out="$2"
            shift
            shift
            ;;

        -p)
            if ! test "$2" -gt 0 2> /dev/null; then
                echo "num-cpu must be positive integer" > /dev/stderr
                exit 1
            fi
            nproc="$2"
            shift
            shift
            ;;

        -r)
            # Set ${run} to a zero-padded, four-digit string
            # representation of the integer.  XXX Rename runno?
            if ! test "$2" -gt 0 2> /dev/null; then
                echo "runno must be positive integer" > /dev/stderr
                exit 1
            fi
            run=`echo "$2" | awk '{ printf("%04d", $1); }'`
            shift
            shift
            ;;

        -x)
            # Experiment subdirectory within the ${CXI}.
            if ! ssh ${NODE} "test -d \"${CXI}/$2\" 2> /dev/null"; then
                echo "exp must be a directory" > /dev/stderr
                exit 1
            fi
            exp="$2"
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
# arguments are present.
if test -z "${cfg}" -o -z "${run}"; then
    echo "Must specify -c and -r options" > /dev/stderr
    exit 1
fi
if test $# -gt 0; then
    echo "Extraneous arguments" > /dev/stderr
    exit 1
fi

# Take ${exp} from the environment unless overridden on the command
# line.
test -n "${EXP}" -a -z "${exp}"&& exp="${EXP}"

# If num-cpu is given in ${cfg}, use that instead of whatever might be
# given on the command line.  Otherwise, the number of processes per
# node should be between 7 and 9 according to Marc Messerschmidt.
t=`awk -F= '/^[[:space:]]*num-cpu[[:space:]]*=/ { \
                printf("%d\n", $2);               \
            }' "${cfg}"`
if test -n "${t}"; then
    test -n "${nproc}" && \
        echo "-p option overridden by configuration file" > /dev/stderr
    nproc="${t}"
elif test -z "${nproc}"; then
    nproc="8"
fi

# Write output to the results subdirectory within the experiment's
# scratch space, unless specified on the command line.  Create the
# directory for the run on the cluster.  Determine the zero-padded,
# three-digit sequence number of the current analysis, by increasing
# the highest existing sequence number by one.
if test -z "${out}"; then
    out="${CXI}/${exp}/scratch/results"
fi
out="${out}/r${run}"
seq=`ssh ${NODE} "mkdir -p \"${out}\" ; ls \"${out}\"" | sort -n | tail -n 1`
if test -z "${seq}" || ! test "${seq}" -ge 0 2> /dev/null; then
    seq="000"
else
    seq=`expr "${seq}" \+ 1`
    seq=`echo "${seq}" | awk '{ printf("%03d", $1); }'`
fi
out="${out}/${seq}"

# Absolute path to the directory with the XTC files.
xtc="${CXI}/${exp}/xtc"

# Sorted list of unique streams for ${run}.
streams=`ssh ${NODE} "ls ${xtc}/e*-r${run}-s*"          \
    | sed -e "s:.*-\(s[[:digit:]][[:digit:]]\)-c.*:\1:" \
    | sort -u                                           \
    | tr -s "\n" " "`

# Create a directory for temporary files, and install a trap to clean
# it all up.
tmpdir=`mktemp -d` || exit 1
trap "ssh ${NODE} \"rm -fr \\\"${out}\\\"\"; \
      rm -fr \"${tmpdir}\";                  \
      exit 1" HUP INT QUIT TERM

# Write a configuration file for the analysis of each stream by
# substituting the directory names with appropriate directories in
# ${out}, and appending the stream number to the base name.  Create a
# script to submit the jobs to the queue.  XXX Dump the environment in
# here, too?
cat > "${tmpdir}/submit.sh" << EOF
#! /bin/sh

NPROC="${nproc}"
OUT="${out}"
PYANA="${PYANA}"
PYTHON="${PYTHON}"
XTC="${xtc}"
EOF
for s in ${streams}; do
    sed -e "s:\([[:alnum:]]\+\)\(_dirname[[:space:]]*=\).*:\1\2 ${out}/\1:"   \
        -e "s:\([[:alnum:]]\+_basename[[:space:]]*=.*\)[[:space:]]*:\1${s}-:" \
        "${cfg}" > "${tmpdir}/pyana_${s}.cfg"

    # Cannot use an indented here-document (<<-), because that would
    # require leading tabs which are not permitted by
    # libtbx.find_clutter.
    cat >> "${tmpdir}/submit.sh" << EOF
bsub -o "\${OUT}/stdout/${s}.out" -q psfehq -J "r${run}-${s}"    \\
    "\"\${PYTHON}\" \"\${PYANA}\" -c \"\${OUT}/pyana_${s}.cfg\" \\
                                -p \"\${NPROC}\"             \\
                                \${XTC}/e*-r${run}-${s}-c*"
EOF
done
chmod 755 "${tmpdir}/submit.sh"

# Create all directories for the output from the analysis.  This
# eliminates a race condition when run in parallel.
directories=`awk -F=                                    \
    '/^[[:space:]]*[[:alnum:]]+_dirname[[:space:]]*=/ { \
         gsub(/^ /, "", $2);                            \
         gsub(/ $/, "", $2);                            \
         printf("\"%s\"\n", $2);                        \
     }' "${tmpdir}"/pyana_s[0-9][0-9].cfg | sort -u | tr -s "\n" " "`
ssh ${NODE} "mkdir -p \"${out}/stdout\" ${directories}"

# Copy the configuration files and the submission script to ${out}.
# Submit the analysis of all streams to the queueing system.
scp -pq "${cfg}"                         "${NODE}:${out}/pyana.cfg"
scp -pq "${tmpdir}"/pyana_s[0-9][0-9].cfg \
        "${tmpdir}/submit.sh"             "${NODE}:${out}"
"${tmpdir}/submit.sh"
rm -fr "${tmpdir}"

echo "Output directory: ${out}"
exit 0
