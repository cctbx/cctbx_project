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

# A random host that has the scratch directory mounted.  psexport is
# preferred over psanafeh, since the latter is not accessible from
# everywhere.
NODE="psexport.slac.stanford.edu"

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

args=`getopt c:o:p:q:r:x: $*`
if test $? -ne 0; then
    echo "Usage: lsf.sh -c config -r runno [-o output] [-p num-cpu] [-q queue] [-x exp]" > /dev/stderr
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

        -q)
            queue="$2"
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

# If no queue is given on the command line then submit to default queue
if [ -z "${queue}" ]; then
    queue="psfehq"
fi

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
# line, and find its absolute path under /reg/data.
test -n "${EXP}" -a -z "${exp}" && exp="${EXP}"
exp=`ssh ${NODE} \
    "find /reg/data -maxdepth 3 -name \"${exp}\" -type d 2> /dev/null"`
if ! test -d "${exp}" 2> /dev/null; then
    echo "Could not find experiment subdirectory for ${exp}" > /dev/stderr
    exit 1
fi

# Construct an absolute path to the directory with the XTC files as
# well as a sorted list of unique stream numbers for ${run}.
xtc="${exp}/xtc"
streams=`ssh ${NODE} "ls ${xtc}/e*-r${run}-s* 2> /dev/null" \
    | sed -e "s:.*-s\([[:digit:]]\+\)-c.*:\1:"              \
    | sort -u                                               \
    | tr -s "\n" " "`
if test -z "${streams}"; then
    echo "No streams in ${xtc}" > /dev/stderr
    exit 1
fi

# If num-cpu is given in ${cfg}, use that instead of whatever might be
# given on the command line.  Otherwise, the number of processes per
# host should be between 7 and 9 according to Marc Messerschmidt.
t=`awk -F= '/^[[:space:]]*num-cpu[[:space:]]*=/ { \
                printf("%d\n", $2);               \
            }' "${cfg}"`
if test "${t}" -gt 0 2> /dev/null; then
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
    out="${exp}/scratch/results"
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

# Create a directory for temporary files, and install a trap to clean
# it all up.
tmpdir=`mktemp -d` || exit 1
trap "ssh ${NODE} \"rm -fr \\\"${out}\\\"\"; \
      rm -fr \"${tmpdir}\";                  \
      exit 1" HUP INT QUIT TERM

# Write a configuration file for the analysis of each stream by
# substituting the directory names with appropriate directories in
# ${out}, and appending the stream number to the base name.  Create a
# run-script for each job, as well as a convenience script to submit
# all the jobs to the queue.  XXX Dump the environment in here, too?
cat > "${tmpdir}/submit.sh" << EOF
#! /bin/sh

OUT="${out}"

EOF
for s in ${streams}; do
    sed -e "s:\([[:alnum:]]\+\)\(_dirname[[:space:]]*=\).*:\1\2 ${out}/\1:"    \
        -e "s:\([[:alnum:]]\+_basename[[:space:]]*=.*\)[[:space:]]*:\1s${s}-:" \
        "${cfg}" > "${tmpdir}/pyana_s${s}.cfg"

    # Process each stream on a single host as a base-1 indexed job.
    # Allocate no more than ${nproc} processors.  Allow the job to
    # start if at least one processor is available on the host.
    # Cannot use an indented here-document (<<-), because that would
    # require leading tabs which are not permitted by
    # libtbx.find_clutter.
    i=`expr "${s}" \+ 1`
    cat >> "${tmpdir}/submit.sh" << EOF
bsub -J "r${run}[${i}]" -n "1,${nproc}" -o "\${OUT}/stdout/s${s}.out" \\
    -q ${queue} -R "span[hosts=1]" "\${OUT}/pyana_s${s}.sh"
EOF
    # limited cores/user:  psfehq.  unlimited: psfehmpiq
    # Create the run-script for stream ${s}.  Fall back on using a
    # single processor if the number of available processors cannot be
    # obtained from the environment.
    cat > "${tmpdir}/pyana_s${s}.sh" << EOF
#! /bin/sh

NPROC=\`printenv LSB_MCPU_HOSTS | awk '{ printf("%d\n", \$2); }'\`

test "\${NPROC}" -gt 0 2> /dev/null || NPROC="1"
"${PYANA}" \\
    -c "${out}/pyana_s${s}.cfg" \\
    -p "\${NPROC}" \\
    "${xtc}"/e*-r${run}-s${s}-c*
EOF
    chmod 755 "${tmpdir}/pyana_s${s}.sh"
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
        "${tmpdir}"/pyana_s[0-9][0-9].sh  \
        "${tmpdir}/submit.sh"             "${NODE}:${out}"
"${tmpdir}/submit.sh"
rm -fr "${tmpdir}"

echo "Output directory: ${out}"
exit 0
