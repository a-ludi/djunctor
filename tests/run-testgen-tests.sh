#!/bin/bash

#-----------------------------------------------------------------------------
# This runs integration tests for the testgen commandline tool.
#
# Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
# License: Subject to the terms of the MIT license, as written in the
#          included LICENSE file.
# Authors: Arne Ludwig <arne.ludwig@posteo.de>
#-----------------------------------------------------------------------------

TEST_DATA_ARCHIVE="testgen.tar.xz"
TEST_DATA_ASSEMBLY1="assembly1"
TEST_DATA_ASSEMBLY2="assembly2"
TEST_DATA_ALIGNMENT="assembly2.assembly1.las"
GDB_INIT_SCRIPT="gdbinit"
TESTGEN_OPTS=(
    -v
    -v
    -v
)
BUILD_OPTS=(--config=testgen)

ARGV=("$@")
KEEP_TEMP=false
RUN_TESTGEN=true
RUN_GDB=false
GDB="${GDB:-gdb}"
if ! [[ -v GDBFLAGS ]]; then
    GDBFLAGS=()
fi
SHOW_COVERAGE=false
SHOW_UNCOVERED_LINES=false
UNCOVERED_LINES_CONTEXT=2
VERBOSE=false
FORCE=false

function init_script()
{
    set -e  # exit on failure
    trap clean_up EXIT

    TEST_ROOT="$(dirname "$(realpath "$0")")"
    LOG_COPY="$PWD/testgen-tests.log"
    RESULTS_ARCHIVE="$PWD/testgen-tests.tar.gz"

    parse_opts

    WORKDIR="$(mktemp --tmpdir -d testgen-tests.XXXXXX)"
    OUTPUT_LOG="$WORKDIR/output.log"
    RESULT_FILE="$WORKDIR/result.fasta"
    RESULT_DB="$WORKDIR/result.dam"
}

function parse_opts()
{
    while getopts "chDfgkuv" OPTION "${ARGV[@]}"; do
        case "$OPTION" in
            c)
                BUILD_OPTS+=('--build=cov')
                SHOW_COVERAGE=true
                ;;
            D)
                RUN_TESTGEN=false
                ;;
            f)
                FORCE=true
                ;;
            g)
                RUN_GDB=true
                ;;
            h)
                usage
                exit
                ;;
            k)
                KEEP_TEMP=true
                TESTGEN_OPTS+=('-k')
                ;;
            u)
                SHOW_UNCOVERED_LINES=true
                ;;
            v)
                VERBOSE=true
                ;;
            *)
                usage
                exit 1
                ;;
        esac
    done
}

function usage()
{
    echo "Usage: ${ARGS[0]} [-cDghu]"
    echo
    echo "Run the integration test suite for testgen."
    echo
    echo "Optional arguments:"
    echo " -c        Enables code coverage statistics to be generated; show coverage"
    echo "           summary after tests."
    echo " -D        Do not run testgen; instead just run tests against the results ($RESULTS_ARCHIVE)."
    echo " -f        Force recreation of generated files"
    echo " -g        Open interactive gdb session and exit afterwards. In gdb type "'`testgen`'
    echo " -h        Prints this help."
    echo " -k        Keep temporary files; this is forwarded to testgen."
    echo " -u[=NUM]  If -c is given report uncovered lines in coverage summary. If given print NUM"
    echo "           lines of context (default: 2)"
    echo " -v        Show more output for debugging."
}

function clean_up()
{
    if $RUN_TESTGEN && ! $RUN_GDB;
    then
        backup_results
        echo "results saved in $RESULTS_ARCHIVE"

        cp "$OUTPUT_LOG" "$LOG_COPY"
        echo "output copied to $LOG_COPY"
    fi

    # remove coverage statistics of libraries
    rm -f -- -*.lst

    if $KEEP_TEMP;
    then
        echo "keeping workdirs; please remove after inspection:"
        echo "    $WORKDIR"
        echo "    $TESTGEN_WORKDIR"
    else
        rm -rf "$WORKDIR"
    fi
}

function backup_results()
{
    local FILE_LIST="$(mktemp --tmpdir testgen-tests.XXXXXX.files.lst)"

    pushd "$WORKDIR" > /dev/null
    find . -type f -print0 > "$FILE_LIST"
    tar -czf "$RESULTS_ARCHIVE" \
        --verbatim-files-from \
        --null \
        -T "$FILE_LIST"
    popd > /dev/null

    rm -f "$FILE_LIST"
}

function restore_results()
{
    pushd "$WORKDIR" > /dev/null
    tar -xf "$RESULTS_ARCHIVE"
    popd > /dev/null
    set_test_data_paths
}

function provide_test_data()
{
    pushd "$WORKDIR" > /dev/null
    cp "$TEST_ROOT/data/$TEST_DATA_ARCHIVE" ./
    tar -xf "$TEST_DATA_ARCHIVE"
    rm "$TEST_DATA_ARCHIVE"
    popd > /dev/null
    set_test_data_paths
}

function set_test_data_paths()
{
    TEST_DATA_ASSEMBLY1="$WORKDIR/$TEST_DATA_ASSEMBLY1"
    TEST_DATA_ASSEMBLY2="$WORKDIR/$TEST_DATA_ASSEMBLY2"
    TEST_DATA_ALIGNMENT="$WORKDIR/$TEST_DATA_ALIGNMENT"
}

function run_testgen()
{
    TESTGEN_OPTS+=("translocate")
    TESTGEN_OPTS+=("$TEST_DATA_ASSEMBLY1.dam")
    TESTGEN_OPTS+=("$TEST_DATA_ASSEMBLY2.dam")
    TESTGEN_OPTS+=("$TEST_DATA_ALIGNMENT")

    if $RUN_GDB;
    then
        build_gdb_init_script > "$WORKDIR/$GDB_INIT_SCRIPT"

        "$GDB" "${GDBFLAGS[@]}" -x "$WORKDIR/$GDB_INIT_SCRIPT" testgen
        exit
    else
        ./testgen "${TESTGEN_OPTS[@]}" \
                        2> "$OUTPUT_LOG" \
                        1> "$RESULT_FILE"
        local TESTGEN_STATUS="$?"
        TESTGEN_WORKDIR="$(json_log | jq -r 'select(has("workdir")) | .workdir')"

        if (( TESTGEN_STATUS != 0 ));
        then
            echo "Error while executing testgen ($TESTGEN_STATUS); see log for details." >&2

            exit 1
        fi
    fi
}

function build_gdb_init_script()
{
    echo 'define testgen'
    echo run "${TESTGEN_OPTS[@]}"
    echo 'end'
    echo
    echo 'echo ------------------------------------\n'
    echo 'echo type `testgen` to start the program\n'
    echo 'echo ------------------------------------\n'
}

function prepare_tests()
{
    pushd "$WORKDIR" > /dev/null
    DAM2fasta -w50 "$TEST_DATA_ASSEMBLY1.dam"
    popd > /dev/null
}

function do_tests()
{
    local NUM_TEST_CASES="$(list_test_cases | wc -l)"
    local TEST_CASE_LOG="$WORKDIR/test-case.log"
    local FAILURES=()

    echo -n "Running $NUM_TEST_CASES test cases: "

    for test_case in $(list_test_cases);
    do
        if $test_case > "$TEST_CASE_LOG";
        then
            if $VERBOSE;
            then
                echo "$test_case: passed"
            else
                echo -n "."
            fi
        else
            if $VERBOSE;
            then
                echo "$test_case: FAILED"
            else
                echo -n "f"
            fi
            FAILURES+=("$test_case: $(cat "$TEST_CASE_LOG")")
        fi
    done

    echo  # finish line

    if (( ${#FAILURES[*]} > 0 ));
    then
        echo "${#FAILURES[*]} cases failed:"

        for failed_test_case in "${FAILURES[@]}";
        do
            echo "    $failed_test_case"
        done

    fi

    echo
    echo "successes: $(( $NUM_TEST_CASES - ${#FAILURES[*]} )) failures: ${#FAILURES[*]} total: $NUM_TEST_CASES"
}

function list_test_cases()
{
    declare -F | grep -oP "(?<=^declare -f )test_.*$" | sort
}

function show_run_time_summary()
{
    local RUN_TIME="$(json_log 'function' | jq -rs 'map(select(.function == "run") | .timestamp) | (.[1] - .[0])/10000000')"

    echo "testgen run time: $RUN_TIME seconds"
}

function show_coverage_summary()
{
    echo "Coverage percentages per file:"
    grep -hoP '^.*is \d+% covered$' *.lst | indent

    if $SHOW_UNCOVERED_LINES;
    then
        echo
        echo "Uncovered lines:"
        grep --context=$UNCOVERED_LINES_CONTEXT -P '^0000000\|' *.lst | indent
        echo $UNCOVERED_LINES_CONTEXT
    fi
}

function indent()
{
    sed 's/^/  /'
}

function main()
{
    init_script
    if $RUN_TESTGEN;
    then
        dub build "${BUILD_OPTS[@]}"
        provide_test_data
        run_testgen

        prepare_tests
    else
        restore_results
    fi
    show_run_time_summary
    do_tests

    if $SHOW_COVERAGE;
    then
        show_coverage_summary
    fi
}


#-----------------------------------------------------------------------------
# Test Helpers
#-----------------------------------------------------------------------------

function expect_json()
{
    local OBSERVED="$(json_log "$4" | jq --sort-keys --slurp "$JQ_DEFS map(select($1))")"

    if ! jq --exit-status "$JQ_DEFS $2" > /dev/null <<<"$OBSERVED";
    then
        echo "expected: $2"
        if $VERBOSE; then
            echo "observed: $OBSERVED"
        fi

        if [[ -n "$3" ]];
        then
            echo "debug: $(jq "$3" <<< "$OBSERVED")"
        fi

        return 1
    fi
}

function json_log()
{
    if [[ -z "$1" ]]; then
        grep --text '^{' "$OUTPUT_LOG"
    else
        grep --text '^{' "$OUTPUT_LOG" | grep "$1"
    fi

}


#-----------------------------------------------------------------------------
# Test Cases
#-----------------------------------------------------------------------------

function test_output_has_correct_gaps()
{
    grep -v '^>' "$RESULT_FILE" | md5sum -  > "$RESULT_FILE.md5"
    if ! (grep -v '^>' "$TEST_DATA_ASSEMBLY1.fasta" | md5sum -c "$RESULT_FILE.md5");
    then
        echo "unexpected result"

        return 1
    fi
}


#-----------------------------------------------------------------------------
# run main (keep at end of file)
#-----------------------------------------------------------------------------

main
