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
GDB_INIT_SCRIPT="gdbinit"
DJUNCTOR_OPTS=(-v -v -v --input-provide-method symlink)
BUILD_OPTS=(--build=debug)
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
    pushd "$WORKDIR" > /dev/null
    cp "$TEST_ROOT/data/$TEST_DATA_ARCHIVE" ./
    tar -xzf "$TEST_DATA_ARCHIVE"
    rm "$TEST_DATA_ARCHIVE"
    popd > /dev/null

    TEST_DATA_READS="$WORKDIR/$TEST_DATA_READS"
    TEST_DATA_REF="$WORKDIR/$TEST_DATA_REF"
    TEST_DATA_MODREF="$WORKDIR/$TEST_DATA_MODREF"
}

function run_djunctor()
{
    if $RUN_GDB;
    then
        build_gdb_init_script > "$WORKDIR/$GDB_INIT_SCRIPT"
        echo '------------------------------------'
        echo 'type `djunctor` to start the program'
        echo '------------------------------------'

        gdb -x "$WORKDIR/$GDB_INIT_SCRIPT" djunctor
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

function build_gdb_init_script()
{
    echo 'define djunctor'
    echo run "${DJUNCTOR_OPTS[@]}" "$TEST_DATA_MODREF.dam" "$TEST_DATA_READS.dam"
    echo 'end'
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
    declare -F | grep -oP "(?<=^declare -f )test_.*$" | sort
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

function expect_insert_sequence_ends()
{
    local BEGIN_CONTIG="$1"
    local END_CONTIG="$2"
    local EXPECTED_HEAD="$3"
    local EXPECTED_TAIL="$4"

    local INSERTION_SEQUENCE="$(json_log | jq --raw-output 'select(has("insertSequence") and .contigIds == ['"$BEGIN_CONTIG"', '"$END_CONTIG"']) | .insertSequence')"
    local INSERTION_HEAD="$(head -c${#EXPECTED_HEAD} <<<"$INSERTION_SEQUENCE")"
    local INSERTION_TAIL="$(head -c-1 <<<"$INSERTION_SEQUENCE" | tail -c${#EXPECTED_TAIL})"

    if [[ "$INSERTION_HEAD" != "$EXPECTED_HEAD" ]];
    then
        echo "insertion head mismatch:"
        echo "  expected: $EXPECTED_HEAD"
        echo "       got: $INSERTION_HEAD"

        return 1
    fi

    if [[ "$INSERTION_TAIL" != "$EXPECTED_TAIL" ]];
    then
        echo "insertion tail mismatch:"
        echo "  expected: $EXPECTED_TAIL"
        echo "       got: $INSERTION_TAIL"

        return 1
    fi
}

#-----------------------------------------------------------------------------
# Test Cases
#-----------------------------------------------------------------------------

function test_front_extension_found()
{
    expect_json \
        '. | has("gapInfo") and .step == "findHits" and .readState == "raw"' \
        '.[0].gapInfo | map(select(.type == "front")) | length == 1 and map(.contigIds) == [[5]] and map(.estimateLengthMean) == [4538] and map(.numReads) == [15]' \
        '.[0].gapInfo | map(select(.type == "front"))' && \
    expect_json \
        '. | has("gapInfo") and .step == "findHits" and .readState == "consensus"' \
        '(map(select(.gapInfo | length == 1 and .[0].type == "front") | .gapInfo) | flatten) | length == 1 and map(.contigIds) == [[5]] and map(.estimateLengthMean) == [4217] and map(.numReads) == [15]' \
        '(map(select(.gapInfo | length == 1 and .[0].type == "front") | .gapInfo) | flatten)'
}

function test_gaps_found()
{
    expect_json \
        '. | has("gapInfo") and .step == "findHits" and .readState == "raw"' \
        '.[0].gapInfo | map(select(.type == "gap")) | length == 3 and map(.contigIds) == [[1, 2], [2, 3], [3, 4]] and map(.estimateLengthMean) == [4156, 8561, 5064] and map(.numReads) == [71, 81, 64]' \
        '.[0].gapInfo | map(select(.type == "gap"))' && \
    expect_json \
        '. | has("gapInfo") and .step == "findHits" and .readState == "consensus"' \
        '(map(select(.gapInfo | length == 1 and .[0].type == "gap") | .gapInfo) | flatten) | length == 3 and map(.contigIds) == [[1, 2], [2, 3], [3, 4]] and map(.estimateLengthMean) == [4100, 8449, 5000] and map(.numReads) == [71, 81, 64]' \
        '(map(select(.gapInfo | length == 1 and .[0].type == "gap") | .gapInfo) | flatten)'
}

function test_back_extension_found()
{
    expect_json \
        '. | has("gapInfo") and .step == "findHits" and .readState == "raw"' \
        '.[0].gapInfo | map(select(.type == "back")) | length == 0' \
        '.[0].gapInfo | map(select(.type == "back"))' && \
    expect_json \
        '. | has("gapInfo") and .step == "findHits" and .readState == "consensus"' \
        '(map(select(.gapInfo | length == 1 and .[0].type == "back") | .gapInfo) | flatten) | length == 0' \
        '(map(select(.gapInfo | length == 1 and .[0].type == "back") | .gapInfo) | flatten)'
}

function test_gaps_filled()
{
    expect_json \
        '.type == "gap" and .step == "insertHits" and has("readId")' \
        'length == 3 and map(.readId) == [17, 6, 2] and map(.contigIds) == [[1, 2], [2, 3], [3, 4]]' \
        '{ readIds: map(.readId), contigIds: map(.contigIds) }'
}

function test_gap_1_2_nicely_filled()
{
    expect_json \
        '.contigIds == [1, 2] and has("transformedInsertionInfo")' \
        '.[0] | ((.transformedInsertionInfo | (.begin.idx == 8300 and .end.idx == 0)) and .transformedInsertionInfo.length == (.insertSequence | length))' \
        '.[0] | {".begin.idx": .transformedInsertionInfo.begin.idx, ".end.idx": .transformedInsertionInfo.end.idx, ".length": .transformedInsertionInfo.length, "length": (.insertSequence | length)}' && \
    expect_insert_sequence_ends \
        1 2 \
        'ctgcagcagagaccacacattaactcttattacgttccaacaacctataa' \
        'acttatgtttactagttattatggtctcaatgtgtgtacacccccacccc'
    # this is the expected consensus
}

function test_gap_2_3_nicely_filled()
{
    expect_json \
        '.contigIds == [2, 3] and has("transformedInsertionInfo")' \
        '.[0].transformedInsertionInfo | (.begin.idx == 20750 and .end.idx == 0 and .length == 8450)' \
        '.[0].transformedInsertionInfo | {".begin.idx": .begin.idx, ".end.idx": .end.idx, ".length": .length}' && \
    expect_insert_sequence_ends \
        2 3 \
        'actccgcttaaatttattttacctaaaatccttttaacagacaaagtcca' \
        'ctaaaagtatactggcttaaatgctaatttcatctgaaaaataccacaga'
    # this is the expected consensus
}

function test_gap_3_4_nicely_filled()
{
    expect_json \
        '.contigIds == [3, 4] and has("transformedInsertionInfo")' \
        '.[0].transformedInsertionInfo | (.begin.idx == 154900 and .end.idx == 0 and .length == 5000)' \
        '.[0].transformedInsertionInfo | {".begin.idx": .begin.idx, ".end.idx": .end.idx, ".length": .length}' && \
    expect_insert_sequence_ends \
        3 4 \
        'cctcaggccggtttgaagaccagcacagcccatgtgagcagggcacaggg' \
        'ttctccagcatcatcttgggagaatgttgacgtcactgatcttggtttct'
    # this is the true sequence
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
        1 $((42 - 0 + 8450 + 8350 - 0 + 4100 + 8300)) && \
    expect_transformed_coord "$COORD_TRANSFORM" \
        4 42 \
        1 $((42 - 0 + 5000 + 125700 - 0 + 8450 + 8350 - 0 + 4100 + 8300)) && \
    expect_transformed_coord "$COORD_TRANSFORM" \
        1337 42 \
        1337 42
}



#-----------------------------------------------------------------------------
# run main (keep at end of file)
#-----------------------------------------------------------------------------

main
