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
GDB_INIT_SCRIPT="gdbinit"
DJUNCTOR_OPTS=(
    -v
    -v
    -v
    --input-provide-method symlink
    --repeat-mask djunctor_repeat
    --reads-error 0.05
    --reference-error 0.01
)
BUILD_OPTS=(--build=debug)
RESULT_TO_REFERENCE_DALIGNER_OPTS=(
    -e.95
    -l4000
    -h128
)
JQ_DEFS='
    def abs: if . < 0 then -1 * . else . end;
    def zip(arrB):
        . as $arrA |
        [range(0; ([($arrA | length), (arrB | length)] | min))] |
        map(
            . as $i |
            [$arrA[$i], arrB[$i]]
        );
    def gapFilled(expectation):
        .beginContig == expectation.beginContig and
        .endContig == expectation.endContig and
        (.gapEnd - expectation.gapEnd | abs) <= 16;
'

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
FORCE=false

function init_script()
{
    set -e  # exit on failure
    trap clean_up EXIT

    TEST_ROOT="$(dirname "$(realpath "$0")")"
    LOG_COPY="$PWD/integration-tests.log"
    RESULTS_ARCHIVE="$PWD/integration-tests.tar.gz"

    parse_opts

    WORKDIR="$(mktemp --tmpdir -d djunctor-integration-tests.XXXXXX)"
    OUTPUT_LOG="$WORKDIR/output.log"
    RESULT_FILE="$WORKDIR/result.fasta"
    RESULT_DB="$WORKDIR/result.dam"
}

function parse_opts()
{
    while getopts "chDfgkuv" OPTION "${ARGV[@]}"; do
        case "$OPTION" in
            c)
                BUILD_OPTS[${#BUILD_OPTS[*]}]='--build=cov'
                SHOW_COVERAGE=true
                ;;
            D)
                RUN_DJUNCTOR=false
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
    echo " -f        Force recreation of generated files"
    echo " -g        Open interactive gdb session and exit afterwards. In gdb type "'`djunctor`'
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
    find . -type f -print0 > "$FILE_LIST"
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
    daligner "${RESULT_TO_REFERENCE_DALIGNER_OPTS}" "$TEST_DATA_REF.dam" "$RESULT_DB"
    daligner "${RESULT_TO_REFERENCE_DALIGNER_OPTS}" "$TEST_DATA_MODREF.dam" "$RESULT_DB"
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

        if $VERBOSE;
        then
            reference_to_result_alignments | \
                awk -F ',' '{
                    if ($1 == 1 && $2 == '"$RESULT_CONTIG"')
                    {
                        printf "        → (%d - (%d - %d)) == %d < '"$MAX_LENGTH_DIFF"' && %d < '"$MAX_NUM_DIFFS"')\n", $6, $10, $9, ($6 - ($10 - $9)), $11
                    }
                }'
        fi

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

function join()
{
    local IFS_BACKUP="$IFS"

    IFS="$1"
    shift
    echo "$*"

    IFS="$IFS_BACKUP"
}


function align_classified_reads_against_reference_mod()
{
    local CLASS="$1"
    local READS_PATH="$WORKDIR/$CLASS"

    if [[ -f "$READS_PATH.dam" ]] && $FORCE; then
        DBrm "$READS_PATH.dam"
    fi

    if [[ ! -f "$READS_PATH.dam" ]];
    then
        local READ_IDS=($(jq '. | map(.'"$CLASS"') | flatten | unique | sort | .[]'))
        local DAMAPPER_CMD=(
            $(json_log 'damapper' | jq --slurp --raw-output 'map(select(has("command") and .command[0] == "damapper")) | .[0].command[0:-2][]')
            "$TEST_DATA_MODREF.dam"
            "$READS_PATH.dam"
        )

        pushd "$WORKDIR" > /dev/null
        DBshow "$TEST_DATA_READS.dam" "${READ_IDS[@]}" | \
            tee "$READS_PATH.fasta" | \
            fasta2DAM -i "$READS_PATH.dam" && \
            "${DAMAPPER_CMD[@]}"
        popd > /dev/null
    fi
}

#-----------------------------------------------------------------------------
# Test Cases
#-----------------------------------------------------------------------------

function test_short_extension_skipped()
{
    expect_json \
        '. | has("extensionLength") and .step == "findHits" and .readState == "raw"' \
        '. | map(
            .effectedContig | inside([1]) and
            .extensionLength < 100
        ) | all' \
        '. | map({effectedContig, extensionLength, extensionType})' \
        'extensionLength'
}

FRONT_EXTENSION_8_LENGTH=2953

function test_front_extension_found()
{
    expect_json \
        '. | has("gapInfo") and .step == "findHits" and .readState == "raw"' \
        '.[0].gapInfo | map(select(.type == "front")) | map(.contigIds) == [
            [1],
            [8]
        ]' \
        '.[0].gapInfo | map(select(.type == "front") | .contigIds)' \
        'gapInfo'
}

function test_front_extensions_filled()
{
    local INSINFO=(
        '{beginContig: 8, endContig: 8, extensionLength: '"$FRONT_EXTENSION_8_LENGTH"'}'
    )

    expect_json \
        '.type == "front" and .step == "insertHits" and has("insertionInfo")' \
        '
            map({
                beginContig: .insertionInfo.begin.contigId,
                endContig: .insertionInfo.end.contigId,
                extensionLength: (.insertionInfo.length - .insertionInfo.end.idx),
            }) |
            zip(['"$(join ',' "${INSINFO[@]}")"']) |
            map(
                .[0] as $got |
                .[1] as $exp |
                $got == $exp
            ) |
            all
        ' \
        '{ insertionInfos: map({ beginContig: .insertionInfo.begin.contigId, endContig: .insertionInfo.end.contigId, extensionLength: (.insertionInfo.length - .insertionInfo.end.idx)}) }' \
        'insertionInfo'
}

function test_back_extension_found()
{
    expect_json \
        '. | has("gapInfo") and .step == "findHits" and .readState == "raw"' \
        '.[0].gapInfo | map(select(.type == "back")) | map(.contigIds) == [
            [7]
        ]' \
        '.[0].gapInfo | map(select(.type == "back") | .contigIds)' \
        'gapInfo'
}

function test_back_extensions_filled()
{
    local INSINFO=(
        '{beginContig: 7, endContig: 7, extensionLength: 11842}'
    )

    expect_json \
        '.type == "back" and .step == "insertHits" and has("insertionInfo")' \
        '
            map({
                beginContig: .insertionInfo.begin.contigId,
                endContig: .insertionInfo.end.contigId,
                extensionLength: .insertionInfo.length,
            }) |
            zip(['"$(join ',' "${INSINFO[@]}")"']) |
            map(
                .[0] as $got |
                .[1] as $exp |
                $got == $exp
            ) |
            all
        ' \
        '{ insertionInfos: map({ beginContig: .insertionInfo.begin.contigId, endContig: .insertionInfo.end.contigId, extensionLength: .insertionInfo.length}) }' \
        'insertionInfo'
}

function test_gaps_found()
{
    expect_json \
        '. | has("gapInfo") and .step == "findHits" and .readState == "raw"' \
        '.[0].gapInfo | map(select(.type == "gap") | .contigIds) == [
            [ 1,  2],
            [ 2,  3],
            [ 3,  4],
            [ 4,  5],
            [ 5,  6],
            [ 6,  7],
            [ 8,  9],
            [ 9, 10],
            [10, 11]
        ]' \
        '.[0].gapInfo | map(select(.type == "gap") | .contigIds)' \
        'gapInfo'
}

function test_gaps_filled()
{
    local INSINFO=(
        '{beginContig: 1, endContig: 2, gapEnd:  12401}'  #   8300 + 4100 =  12400
        '{beginContig: 2, endContig: 3, gapEnd:  16800}'  #   8350 + 8450 =  16800
        '{beginContig: 3, endContig: 4, gapEnd: 130701}'  # 125700 + 5000 = 130700
        '{beginContig: 4, endContig: 5, gapEnd:  15006}'  #  10000 + 5000 =  15000
        '{beginContig: 5, endContig: 6, gapEnd:  28750}'  #  25750 + 3000 =  28750
        '{beginContig: 6, endContig: 7, gapEnd:  15250}'  #  12750 + 2500 =  15250
        '{beginContig: 8, endContig: 9, gapEnd:  24900}'  #  21400 + 3500 =  24900
    )

    expect_json \
        '.type == "gap" and .step == "insertHits" and has("readId")' \
        '
            map({
                beginContig: .insertionInfo.begin.contigId,
                endContig: .insertionInfo.end.contigId,
                gapEnd: (.insertionInfo.begin.idx + .insertionInfo.length - .insertionInfo.end.idx),
            }) |
            zip(['"$(join ',' "${INSINFO[@]}")"']) |
            map(
                .[0] as $got |
                .[1] as $exp |
                $got | gapFilled($exp)
            ) |
            all
        ' \
        '{ insertionInfos: map({beginContig: .insertionInfo.begin.contigId, endContig: .insertionInfo.end.contigId, gapEnd: (.insertionInfo.begin.idx + .insertionInfo.length - .insertionInfo.end.idx)}) }' \
        'insertionInfo'
}


function test_pile_ups_contain_enough_valid_read()
{
    local MIN_CORRECT_READS_RATIO="0.9"
    local COMPUTED_PILE_UPS="$(json_log 'pileUps' | \
        jq -c 'select(has("pileUps")) | .pileUps | map({ type: .type, contigs: (.readAlignments | map(.[0].contigB.id) | unique | sort), reads: (.readAlignments | map(.[0].contigA.id) | unique | sort) })')"
    local TRUE_PILE_UPS="$(jq -c 'map({ contigs, reads: (.reads | map(.readId) | sort | unique) })' < "$WORKDIR/pile_ups.json")"
    local TEST_RESULT="$(jq '
        def correctReadsRatio: (.correctReads | length) / (.computedReads | length);
        def trueReadsFoundRatio: (.correctReads | length) / (.trueReads | length);
        def minCorrectReadsRatio: '"$MIN_CORRECT_READS_RATIO"';
        def isGoodPileUp: correctReadsRatio >= minCorrectReadsRatio;
        def isBadPileUp: isGoodPileUp | not;
        def buildOutput: {
            contigs,
            type,
            isGoodEnough: isGoodPileUp,
            correctReads: .correctReads,
            incorrectReads: (.computedReads - .correctReads),
            unusedReads: (.trueReads - .computedReads),
            correctReadsRatio: correctReadsRatio,
            trueReadsFoundRatio: trueReadsFoundRatio,
        };

        .truePileUps as $truePileUps |
        .computedPileUps as $computedPileUps |
        $computedPileUps |
        map(
            .reads as $computedReads |
            .type as $type |
            .contigs as $contigs |
            $truePileUps |
            map(
                .reads as $trueReads |
                select(
                    ($type == "front" and .contigs[1] == $contigs[0]) or
                    ($type == "gap" and .contigs == $contigs) or
                    ($type == "back" and .contigs[0] == $contigs[0])
                ) |
                ($computedReads - ($computedReads - $trueReads)) as $correctReads |
                {
                    contigs: $contigs,
                    type: $type,
                    correctReads: $correctReads,
                    computedReads: $computedReads,
                    trueReads: $trueReads,
                }
            ) |
            select(length > 0)
        ) |
        flatten |
        map(buildOutput)
    ' <<<"{\"truePileUps\":$TRUE_PILE_UPS,\"computedPileUps\":$COMPUTED_PILE_UPS}")"

    if $VERBOSE; then
        align_classified_reads_against_reference_mod 'incorrectReads' <<<"$TEST_RESULT"
        align_classified_reads_against_reference_mod 'unusedReads' <<<"$TEST_RESULT"
    fi

    if ! jq --exit-status 'map(select(.isGoodEnough | not)) | length == 0' <<<"$TEST_RESULT" &> /dev/null; then
        local NUM_PILE_UPS="$(jq -r 'length' <<<"$TEST_RESULT")"
        local NUM_GOOD_PILE_UPS="$(jq -r 'map(select(.isGoodEnough)) | length' <<<"$TEST_RESULT")"

        echo -n "found $NUM_GOOD_PILE_UPS/$NUM_PILE_UPS good pile ups"

        if $VERBOSE; then
            echo -n ": "
            jq '.' <<<"$TEST_RESULT"
        else
            echo
        fi

        return 1
    fi
}

function test_result_contigs_properly_align_to_reference()
{
    local NUM_RESULT_CONTIGS="$(json_log 'numReferenceContigs' | \
        jq --raw-output 'select(has("numReferenceContigs")) | .contigSources | length')"

    for (( I = 1; I <= NUM_RESULT_CONTIGS; I++ ));
    do
        result_contig_properly_aligns_to_reference "$I" || return 1
    done
}

function test_number_of_iterations()
{
    expect_json \
        '. | has("iteration")' \
        '. | max_by(.iteration).iteration == 0' \
        'map(.iteration)' \
        'iteration'
}


#-----------------------------------------------------------------------------
# run main (keep at end of file)
#-----------------------------------------------------------------------------

main
