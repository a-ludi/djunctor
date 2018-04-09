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
DJUNCTOR_OPTS=(-v -v -v)

function init_script()
{
    # set -x  # echo commands
    set -e  # exit on failure
    trap clean_up EXIT

    TEST_ROOT="$(dirname "$(realpath "$0")")"
    WORKDIR="$(mktemp --tmpdir -d djunctor-integration-tests.XXXXXX)"
    OUTPUT_LOG="$WORKDIR/output.log"
}

function clean_up()
{
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
    ./djunctor "${DJUNCTOR_OPTS[@]}" \
               "$TEST_DATA_MODREF.dam" \
               "$TEST_DATA_READS.dam" \
               &> "$OUTPUT_LOG"
}

function do_tests()
{
    local NUM_TEST_CASES="$(list_test_cases | wc -l)"
    local TEST_CASE_LOG="$WORKDIR/test-case.log"
    local FAILURES=()

    echo -n "Running $NUM_TEST_CASES test cases: "

    for test_case in $(list_test_cases);
    do
        if $test_case &> "$TEST_CASE_LOG";
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

        cp "$OUTPUT_LOG" ./integration-tests.log
        echo
        echo "output copied to ./integration-tests.log"
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

function expect_line_present()
{
    if ! grep -qP "$1" < "$OUTPUT_LOG"; then
        echo "expected '$1' to be in the output"

        false
    fi
}


#-----------------------------------------------------------------------------
# Test Cases
#-----------------------------------------------------------------------------

function test_extension_found()
{
    expect_line_present 'extending 5 contigs'
    expect_line_present 'extending size \(raw averages\): \[2639, 4455, 3611, 5385, 1306\]'
}

function test_gaps_found()
{
    expect_line_present 'spanning 3 gaps'
    expect_line_present 'spanning size \(raw averages\): \[4154, 8563, 5082\]'
}

function test_number_of_iterations()
{
    expect_line_present 'i: 0'
    ! expect_line_present 'i: 1'
}


#-----------------------------------------------------------------------------
# run main (keep at end of file)
#-----------------------------------------------------------------------------

main
