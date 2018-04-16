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
DJUNCTOR_OPTS=(-v -v -v --input-provide-method copy)
JQ_DEFS=''

function init_script()
{
    # set -x  # echo commands
    set -e  # exit on failure
    trap clean_up EXIT

    TEST_ROOT="$(dirname "$(realpath "$0")")"
    WORKDIR="$(mktemp --tmpdir -d djunctor-integration-tests.XXXXXX)"
    OUTPUT_LOG="$WORKDIR/output.log"
    LOG_COPY="./integration-tests.log"
}

function clean_up()
{
    if [[ -f "$OUTPUT_LOG" ]];
    then
        cp "$OUTPUT_LOG" "$LOG_COPY"
        echo
        echo "output copied to $LOG_COPY"
    fi

    rm -rf "$WORKDIR"
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
    if ! ./djunctor "${DJUNCTOR_OPTS[@]}" \
                    "$TEST_DATA_MODREF.dam" \
                    "$TEST_DATA_READS.dam" \
                    &> "$OUTPUT_LOG";
    then
        echo "Error while executing djunctor ($?); see log for details." >&2

        exit 1
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

function main()
{
    init_script
    dub build
    provide_test_data
    run_djunctor
    do_tests
}


#-----------------------------------------------------------------------------
# Test Helpers
#-----------------------------------------------------------------------------

function expect_json()
{
    local OBSERVED="$(grep '^{' "$OUTPUT_LOG" | jq --sort-keys --slurp "$JQ_DEFS map(select($1))")"

    if ! jq --exit-status "$JQ_DEFS $2" &> /dev/null <<<"$OBSERVED";
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

function test_number_of_iterations()
{
    expect_json \
        '. | has("iteration")' \
        '. | max_by(.iteration).iteration == 0'
}


#-----------------------------------------------------------------------------
# run main (keep at end of file)
#-----------------------------------------------------------------------------

main
