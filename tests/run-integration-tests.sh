#!/bin/bash

#-----------------------------------------------------------------------------
# This runs integration tests for the djuntor commandline tool.
#
# Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
# License: Subject to the terms of the MIT license, as written in the
#          included LICENSE file.
# Authors: Arne Ludwig <arne.ludwig@posteo.de>
#-----------------------------------------------------------------------------

TEST_DATA_ARCHIVE="integration-tests.tar.gz"
TEST_DATA_READS="integration-tests/reads"
TEST_DATA_REF="integration-tests/reference"
TEST_DATA_MODREF="integration-tests/reference_mod"
DJUNCTOR_OPTS=(-v -v -v --input-provide-method symlink)
BUILD_OPTS=()
JQ_DEFS=''

ARGV=("$@")
KEEP_TEMP=false
RUN_DJUNCTOR=true
RUN_GDB=false
SHOW_COVERAGE=false
SHOW_UNCOVERED_LINES=false
UNCOVERED_LINES_CONTEXT=2

function init_script()
{
    set -e  # exit on failure
    trap clean_up EXIT

    TEST_ROOT="$(dirname "$(realpath "$0")")"
    LOG_COPY="./integration-tests.log"

    parse_opts

    WORKDIR="$(mktemp --tmpdir -d djunctor-integration-tests.XXXXXX)"

    if $RUN_DJUNCTOR;
    then
        OUTPUT_LOG="$WORKDIR/output.log"
    else
        OUTPUT_LOG="$LOG_COPY"
    fi
}

function parse_opts()
{
    while getopts "chDgTu" OPTION "${ARGV[@]}"; do
        case "$OPTION" in
            c)
                BUILD_OPTS[${#BUILD_OPTS[*]}]='--build=cov'
                SHOW_COVERAGE=true
                ;;
            D)
                RUN_DJUNCTOR=false
                ;;
            g)
                RUN_GDB=true
                ;;
            h)
                usage
                exit
                ;;
            T)
                KEEP_TEMP=true
                DJUNCTOR_OPTS[${#DJUNCTOR_OPTS[*]}]='-T'
                ;;
            u)
                SHOW_UNCOVERED_LINES=true
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
    echo "Run the integration test suite for djunctor."
    echo
    echo "Optional arguments:"
    echo " -c        Enables code coverage statistics to be generated; show coverage"
    echo "           summary after tests."
    echo " -D        Do not run djunctor; instead just run tests against the last log ($LOG_COPY)."
    echo " -g        Open interactive gdb session and exit afterwards. Prints the "'`run`'" command"
    echo "           to be used in gdb"
    echo " -h        Prints this help."
    echo " -T        Keep temporary files; this is forwarded to djunctor."
    echo " -u[=NUM]  If -c is given report uncovered lines in coverage summary. If given print NUM"
    echo "           lines of context (default: 2)"
}

function clean_up()
{
    if $RUN_DJUNCTOR && [[ -f "$OUTPUT_LOG" ]];
    then
        cp "$OUTPUT_LOG" "$LOG_COPY"
        echo
        echo "output copied to $LOG_COPY"
    fi

    # remove coverage statistics of libraries
    rm -f -- -*.lst

    if $KEEP_TEMP;
    then
        echo "keeping workdir '$WORKDIR'; please remove after inspection"
    else
        rm -rf "$WORKDIR"
    fi
}

function provide_test_data()
{
    pushd "$WORKDIR"
    cp "$TEST_ROOT/data/$TEST_DATA_ARCHIVE" ./
    tar -xzf "$TEST_DATA_ARCHIVE"
    rm "$TEST_DATA_ARCHIVE"
    popd

    TEST_DATA_READS="$WORKDIR/$TEST_DATA_READS"
    TEST_DATA_REF="$WORKDIR/$TEST_DATA_REF"
    TEST_DATA_MODREF="$WORKDIR/$TEST_DATA_MODREF"
}

function run_djunctor()
{
    if $RUN_GDB;
    then
        echo '------------------------'
        echo run "${DJUNCTOR_OPTS[@]}" \
                   "$TEST_DATA_MODREF.dam" \
                   "$TEST_DATA_READS.dam"
        echo '------------------------'
        gdb djunctor
        exit
    else
        if ! ./djunctor "${DJUNCTOR_OPTS[@]}" \
                        "$TEST_DATA_MODREF.dam" \
                        "$TEST_DATA_READS.dam" \
                        &> "$OUTPUT_LOG";
        then
            echo "Error while executing djunctor ($?); see log for details." >&2

            exit 1
        fi
    fi
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
            echo -n "."
        else
            echo -n "f"
            FAILURES[${#FAILURES[*]}]="$test_case: $(cat "$TEST_CASE_LOG")"
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
}

function list_test_cases()
{
    declare -F | grep -oP "(?<=^declare -f )test_.*$"
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
    if $RUN_DJUNCTOR;
    then
        dub build "${BUILD_OPTS[@]}"
        provide_test_data
        run_djunctor
    fi
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
    local OBSERVED="$(json_log | jq --sort-keys --slurp "$JQ_DEFS map(select($1))")"

    if ! jq --exit-status "$JQ_DEFS $2" > /dev/null <<<"$OBSERVED";
    then
        echo "expected: $2"
        echo "observed: $OBSERVED"

        if [[ -n "$3" ]];
        then
            echo "debug: $(jq "$3" <<< "$OBSERVED")"
        fi

        return 1
    fi
}

function json_log()
{
    grep '^{' "$OUTPUT_LOG"
}

function expect_transformed_coord()
{
    local COORD_TRANSFORM="$1"
    local IN_CONTIG="$2"
    local IN_IDX="$3"
    local OUT_CONTIG="$4"
    local OUT_IDX="$5"

    local OBSERVED="$(python "$COORD_TRANSFORM" "$IN_CONTIG" "$IN_IDX")"
    local EXPECTED="$OUT_CONTIG $OUT_IDX"

    if [[ "$OBSERVED" != "$EXPECTED" ]]; then
        echo "expected: $EXPECTED observed: $OBSERVED"

        return 1
    fi
}

#-----------------------------------------------------------------------------
# Test Cases
#-----------------------------------------------------------------------------

function test_extension_found()
{
    expect_json \
        '. | has("estimateLengths") and .type == "extension" and .step == "findHits" and .readState == "raw"' \
        'length == 1 and .[0].estimateLengths == [1851, 3898, 3249, 5485, 2175] and  .[0].numReads == [20, 55, 54, 42, 28]'
}

function test_gaps_found()
{
    expect_json \
        '. | has("estimateLengths") and .type == "span" and .step == "findHits" and .readState == "raw" and .numGaps > 1' \
        'length == 1 and .[0].estimateLengths == [4156, 8561, 5064] and .[0].numReads == [38, 6, 4]' && \
    expect_json \
        '. | has("estimateLengths") and .type == "span" and .step == "findHits" and .readState == "consensus" and .numGaps == 1' \
        '(map(.estimateLengths) | flatten == [4100, 8450, 5003]) and (map(.numReads) | flatten == [36, 6, 4])'
}

function test_gaps_filled()
{
    expect_json \
        '. | .type == "span" and .step == "insertHits"' \
        'length == 3 and map(.readId) == [1, 1, 1] and map(.contigIds) == [[1, 2], [2, 3], [3, 4]]' \
        '{ readIds: map(.readId), contigIds: map(.contigIds) }'
}

function test_number_of_iterations()
{
    expect_json \
        '. | has("iteration")' \
        '. | max_by(.iteration).iteration == 0'
}

function test_coordinate_transform()
{
    local COORD_TRANSFORM="$WORKDIR/coord_transform.py"

    json_log | jq --raw-output 'select(has("coordTransformPython")) | .coordTransformPython' > "$COORD_TRANSFORM" && \

    expect_transformed_coord "$COORD_TRANSFORM" \
        1 42 \
        1 $((42)) && \
    expect_transformed_coord "$COORD_TRANSFORM" \
        2 42 \
        1 $((42 - 0 + 4100 + 8300)) && \
    expect_transformed_coord "$COORD_TRANSFORM" \
        3 42 \
        1 $((42 - 0 + 8450 + 125700 - 0 + 4100 + 8300)) && \
    expect_transformed_coord "$COORD_TRANSFORM" \
        4 42 \
        1 $((42 - 0 + 5004 + 125700 - 0 + 8450 + 125700 - 0 + 4100 + 8300)) && \
    expect_transformed_coord "$COORD_TRANSFORM" \
        1337 42 \
        1337 42
}



#-----------------------------------------------------------------------------
# run main (keep at end of file)
#-----------------------------------------------------------------------------

main
