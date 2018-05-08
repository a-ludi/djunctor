#!/bin/bash

#-----------------------------------------------------------------------------
# This runs integration tests for the djuntor commandline tool.
#
# Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
# License: Subject to the terms of the MIT license, as written in the
#          included LICENSE file.
# Authors: Arne Ludwig <arne.ludwig@posteo.de>
#-----------------------------------------------------------------------------

TEST_DATA_ARCHIVE="integration-tests.tar.xz"
TEST_DATA_READS="reads"
TEST_DATA_REF="reference"
TEST_DATA_MODREF="reference_mod"
COORD_TRANSFORM="coord-transform.py"
GDB_INIT_SCRIPT="gdbinit"
DJUNCTOR_OPTS=(-v -v -v --input-provide-method symlink)
BUILD_OPTS=(--build=debug)
JQ_DEFS=''

ARGV=("$@")
KEEP_TEMP=false
RUN_DJUNCTOR=true
RUN_GDB=false
GDB="${GDB:-gdb}"
if ! [[ -v GDBFLAGS ]]; then
    GDBFLAGS=()
fi
SHOW_COVERAGE=false
SHOW_UNCOVERED_LINES=false
UNCOVERED_LINES_CONTEXT=2
VERBOSE=false

function init_script()
{
    set -e  # exit on failure
    trap clean_up EXIT

    TEST_ROOT="$(dirname "$(realpath "$0")")"
    LOG_COPY="$PWD/integration-tests.log"
    RESULTS_ARCHIVE="$PWD/integration-tests.tar.gz"

    parse_opts

    WORKDIR="$(mktemp --tmpdir -d djunctor-integration-tests.XXXXXX)"
    DJUNCTOR_OPTS[${#DJUNCTOR_OPTS[*]}]="--coord-transform"
    DJUNCTOR_OPTS[${#DJUNCTOR_OPTS[*]}]="$WORKDIR/$COORD_TRANSFORM"
    OUTPUT_LOG="$WORKDIR/output.log"
    RESULT_FILE="$WORKDIR/result.fasta"
    RESULT_DB="$WORKDIR/result.dam"
    COORD_TRANSFORM="$WORKDIR/$COORD_TRANSFORM"
}

function parse_opts()
{
    while getopts "chDgkuv" OPTION "${ARGV[@]}"; do
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
            k)
                KEEP_TEMP=true
                DJUNCTOR_OPTS[${#DJUNCTOR_OPTS[*]}]='-k'
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
    echo "Run the integration test suite for djunctor."
    echo
    echo "Optional arguments:"
    echo " -c        Enables code coverage statistics to be generated; show coverage"
    echo "           summary after tests."
    echo " -D        Do not run djunctor; instead just run tests against the results ($RESULTS_ARCHIVE)."
    echo " -g        Open interactive gdb session and exit afterwards. Prints the "'`run`'" command"
    echo "           to be used in gdb"
    echo " -h        Prints this help."
    echo " -k        Keep temporary files; this is forwarded to djunctor."
    echo " -u[=NUM]  If -c is given report uncovered lines in coverage summary. If given print NUM"
    echo "           lines of context (default: 2)"
    echo " -v        Show more output for debugging."
}

function clean_up()
{
    if $RUN_DJUNCTOR && ! $RUN_GDB;
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
        echo "    $DJUNCTOR_WORKDIR"
    else
        rm -rf "$WORKDIR"
    fi
}

function backup_results()
{
    local FILE_LIST="$(mktemp --tmpdir djunctor-integration-tests.XXXXXX.files.lst)"

    pushd "$WORKDIR" > /dev/null
    find . -print0 > "$FILE_LIST"
    tar -czf "$RESULTS_ARCHIVE" \
        --exclude='*.fasta' \
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
    TEST_DATA_READS="$WORKDIR/$TEST_DATA_READS"
    TEST_DATA_REF="$WORKDIR/$TEST_DATA_REF"
    TEST_DATA_MODREF="$WORKDIR/$TEST_DATA_MODREF"
}

function run_djunctor()
{
    if $RUN_GDB;
    then
        build_gdb_init_script > "$WORKDIR/$GDB_INIT_SCRIPT"

        "$GDB" "${GDBFLAGS[@]}" -x "$WORKDIR/$GDB_INIT_SCRIPT" djunctor
        exit
    else
        ./djunctor "${DJUNCTOR_OPTS[@]}" \
                        "$TEST_DATA_MODREF.dam" \
                        "$TEST_DATA_READS.dam" \
                        2> "$OUTPUT_LOG" \
                        1> "$RESULT_FILE"
        local DJUNTOR_STATUS="$?"
        DJUNCTOR_WORKDIR="$(json_log | jq -r 'select(has("workdir")) | .workdir')"

        if (( DJUNTOR_STATUS != 0 ));
        then
            echo "Error while executing djunctor ($DJUNTOR_STATUS); see log for details." >&2

            exit 1
        fi
    fi
}

function build_gdb_init_script()
{
    echo 'define djunctor'
    echo run "${DJUNCTOR_OPTS[@]}" "$TEST_DATA_MODREF.dam" "$TEST_DATA_READS.dam"
    echo 'end'
    echo
    echo 'echo ------------------------------------\n'
    echo 'echo type `djunctor` to start the program\n'
    echo 'echo ------------------------------------\n'
}

function prepare_tests()
{
    fasta2DAM "$RESULT_DB" "$RESULT_FILE" && DBsplit "$RESULT_DB"

    pushd "$WORKDIR" > /dev/null
    daligner "$TEST_DATA_REF.dam" "$RESULT_DB"
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

    echo "djunctor run time: $RUN_TIME seconds"
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

function expect_transformed_coord()
{
    local IN_CONTIG="$1"
    local IN_IDX="$2"
    local OUT_CONTIG="$3"
    local OUT_IDX="$4"

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

function result_contig_properly_aligns_to_reference()
{
    local RESULT_CONTIG="$1"
    local MAX_LENGTH_DIFF=16
    local MAX_NUM_DIFFS=256
    local NUM_MATCHING_ALIGNMENTS=$(reference_to_result_alignments | \
        awk -F ',' '
        {
            if ($1 == 1 && $2 == '"$RESULT_CONTIG"' && ($6 - ($10 - $9)) < '"$MAX_LENGTH_DIFF"' && $11 < '"$MAX_NUM_DIFFS"')
            {
                print
            }
        }' | \
        wc -l)
    if ! (( NUM_MATCHING_ALIGNMENTS == 1 )); then
        echo "expected to find proper alignment for contig $RESULT_CONTIG ($NUM_MATCHING_ALIGNMENTS)"
        reference_to_result_alignments | \
            awk -F ',' '{
                if ($1 == 1 && $2 == '"$RESULT_CONTIG"')
                {
                    printf "        → (%d - (%d - %d)) == %d < '"$MAX_LENGTH_DIFF"' && %d < '"$MAX_NUM_DIFFS"')\n", $6, $10, $9, ($6 - ($10 - $9)), $11
                }
            }'

        return 1
    fi
}

function reference_to_result_alignments()
{
    LAdump_csv "$TEST_DATA_REF.dam" "$RESULT_DB" "$WORKDIR/reference.result.las"
}

function LAdump_csv()
{
    LAdump -c -d -l "$@" | \
        tail -n+3 | \
        tr ' \n' ',,' | \
        sed -e 's/,P/\nP/g' -e 's/[PCDL],//g'
}


#-----------------------------------------------------------------------------
# Test Cases
#-----------------------------------------------------------------------------

function test_front_extension_found()
{
    expect_json \
        '. | has("gapInfo") and .step == "findHits" and .readState == "raw"' \
        '.[0].gapInfo | map(select(.type == "front")) | length == 1 and map(.contigIds) == [[5]]' \
        '.[0].gapInfo | map(select(.type == "front"))' \
        'gapInfo'
}

function test_gaps_found()
{
    expect_json \
        '. | has("gapInfo") and .step == "findHits" and .readState == "raw"' \
        '.[0].gapInfo | map(select(.type == "gap")) | length == 3 and map(.contigIds) == [[1, 2], [2, 3], [3, 4]]' \
        '.[0].gapInfo | map(select(.type == "gap"))' \
        'gapInfo'
}

function test_back_extension_found()
{
    expect_json \
        '. | has("gapInfo") and .step == "findHits" and .readState == "raw"' \
        '.[0].gapInfo | map(select(.type == "back")) | length == 0' \
        '.[0].gapInfo | map(select(.type == "back"))' \
        'gapInfo'
}

function test_front_extensions_filled()
{
    local INSINFO1='{begin: {contigId: 5, idx: 0}, end: {contigId: 5, idx: 800}, length: 6013}'

    expect_json \
        '.type == "front" and .step == "insertHits" and has("insertionInfo")' \
        'length == 1 and map(.insertionInfo) == ['"$INSINFO1"']' \
        '{ insertionInfos: map(.insertionInfo) }' \
        'insertionInfo'
}

function test_gaps_filled()
{
    local INSINFO1='{begin: {contigId: 1, idx: 7900}, end: {contigId: 2, idx: 200}, length: 4700}'
    local INSINFO2='{begin: {contigId: 2, idx: 8300}, end: {contigId: 3, idx: 200}, length: 8700}'
    local INSINFO3='{begin: {contigId: 3, idx: 124400}, end: {contigId: 4, idx: 100}, length: 6400}'

    expect_json \
        '.type == "gap" and .step == "insertHits" and has("readId")' \
        'length == 3 and map(.insertionInfo) == ['"$INSINFO1, $INSINFO2, $INSINFO3"']' \
        '{ insertionInfos: map(.insertionInfo) }' \
        'insertionInfo'
}

function test_merged_result()
{
    expect_json \
        '. | has("numReferenceContigs")' \
        '.[0] | .numReferenceContigs == 5 and (.contigSources | map(.[0]) == [1, 5] and .[0][1].dbFile != .[1][1].dbFile)' \
        'map(.numReferenceContigs)' \
        'numReferenceContigs'
}

function test_result_contig_1_properly_aligns_to_reference()
{
    result_contig_properly_aligns_to_reference 1
}

function test_result_contig_2_properly_aligns_to_reference()
{
    result_contig_properly_aligns_to_reference 2
}

function test_number_of_iterations()
{
    expect_json \
        '. | has("iteration")' \
        '. | max_by(.iteration).iteration == 0' \
        'map(.iteration)' \
        'iteration'
}

function test_coordinate_transform__contig_1()
{
    # This is the transformed coordinate after *perfect* gap filling
    expect_transformed_coord \
        1 42 \
        1 $((42))
}

function test_coordinate_transform__contig_2()
{
    # This is the transformed coordinate after *perfect* gap filling
    expect_transformed_coord \
        2 42 \
        1 $((42 + 12400))
}

function test_coordinate_transform__contig_3()
{
    # This is the transformed coordinate after *perfect* gap filling
    expect_transformed_coord \
        3 42 \
        1 $((42 + 29200))
}

function test_coordinate_transform__contig_4()
{
    # This is the transformed coordinate after *perfect* gap filling
    expect_transformed_coord \
        4 42 \
        1 $((42 + 159900))
}

function test_coordinate_transform__contig_5()
{
    expect_transformed_coord \
        5 42 \
        2 $((42 + 6013 - 800))
}



#-----------------------------------------------------------------------------
# run main (keep at end of file)
#-----------------------------------------------------------------------------

main
