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
TEST_DATA_OUT_MASK="djunctor_repeat"
GDB_INIT_SCRIPT="gdbinit"
DJUNCTOR_OPTS=(
    -v
    -v
    -v
    --input-provide-method symlink
    --reads-error 0.05
    --reference-error 0.01
    --join-policy contigs
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
                BUILD_OPTS+=('--build=cov')
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
                DJUNCTOR_OPTS+=('-k')
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
    TEST_DATA_OUT_MASK="$WORKDIR/$TEST_DATA_OUT_MASK"
}

function run_djunctor()
{
    DJUNCTOR_OPTS+=("--out-mask" "$TEST_DATA_OUT_MASK")
    DJUNCTOR_OPTS+=("$TEST_DATA_MODREF.dam")
    DJUNCTOR_OPTS+=("$TEST_DATA_READS.dam")

    if $RUN_GDB;
    then
        build_gdb_init_script > "$WORKDIR/$GDB_INIT_SCRIPT"

        "$GDB" "${GDBFLAGS[@]}" -x "$WORKDIR/$GDB_INIT_SCRIPT" djunctor
        exit
    else
        ./djunctor "${DJUNCTOR_OPTS[@]}" \
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
    echo run "${DJUNCTOR_OPTS[@]}"
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
    # readAId,readBId,complement,chainPart,readALength,readBLength,readABegin,readAEnd,readBBegin,readBEnd,numDiffs
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

function test_small_pile_ups_skipped()
{
    local MIN_READS_PER_PILE_UP="$(json_log 'minReadsPerPileUp' | \
        jq --raw-output 'select(has("minReadsPerPileUp")) | .minReadsPerPileUp')"

    expect_json \
        '. | has("pileUps")' \
        '.[0].pileUps | map(.readAlignments | length >= '"$MIN_READS_PER_PILE_UP"') | all' \
        '.[0].pileUps | map(select(.readAlignments | length < '"$MIN_READS_PER_PILE_UP"'))' \
        'pileUps'
}

function test_insertions_found()
{
    expect_json \
        '. | has("pileUps")' \
        '.[0].pileUps | map({ type: .type, contigIds: (.readAlignments | map(.[0].alignment.contigA.id) | unique | sort) }) == [
            { type: "front", contigIds: [1] },
            { type: "gap",   contigIds: [1, 2] },
            { type: "gap",   contigIds: [2, 3] },
            { type: "gap",   contigIds: [3, 4] },
            { type: "gap",   contigIds: [4, 5] },
            { type: "gap",   contigIds: [5, 6] },
            { type: "gap",   contigIds: [6, 7] },
            { type: "back",  contigIds: [7] },
            { type: "front", contigIds: [8] },
            { type: "back",  contigIds: [8] },
            { type: "front", contigIds: [9] },
            { type: "gap",   contigIds: [9, 10] },
            { type: "gap",   contigIds: [10, 11] },
            { type: "back",  contigIds: [11] },
            { type: "front", contigIds: [12] },
            { type: "gap",   contigIds: [12, 14] },
            { type: "gap",   contigIds: [13, 14] },
            { type: "gap",   contigIds: [13, 15] }
        ]' \
        '.[0].pileUps | map({ type: .type, contigIds: (.readAlignments | map(.[0].alignment.contigA.id) | unique | sort) })' \
        'pileUps'
}

function test_pile_ups_contain_enough_valid_read()
{
    local MIN_CORRECT_READS_RATIO="0.9"
    local COMPUTED_PILE_UPS="$(json_log 'pileUps' | \
        jq -c 'select(has("pileUps")) |
               .pileUps |
               map({
                   type: .type,
                   contigs: (.readAlignments | map(.[0].alignment.contigA.id) | unique | sort),
                   reads: (.readAlignments | map(.[0].alignment.contigB.id) | unique | sort),
               })')"
    local TRUE_PILE_UPS="$(jq -c 'map({ contigs, reads: (.reads | map(.readId) | sort | unique) })' < "$WORKDIR/pile_ups.json")"
    local TEST_RESULT="$(jq '
        def correctReadsRatio: if (.computedReads | length) > 0 then (.correctReads | length) / (.computedReads | length) else 1 end;
        def trueReadsFoundRatio: if (.trueReads | length) > 0 then (.correctReads | length) / (.trueReads | length) else 1 end;
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
        $truePileUps |
        map(
            .reads as $trueReads |
            .contigs as $trueContigs |
            reduce $computedPileUps[] as $computedPileUp (
                {
                    unusedReads: .reads,
                    computedReads: [],
                };
                if (
                    ($computedPileUp.type == "front" and $trueContigs[1] == $computedPileUp.contigs[0]) or
                    ($computedPileUp.type == "gap" and $trueContigs == $computedPileUp.contigs) or
                    ($computedPileUp.type == "back" and $trueContigs[0] == $computedPileUp.contigs[0])
                )
                then
                    {
                        unusedReads: (.unusedReads - $computedPileUp.reads),
                        computedReads: (.computedReads + $computedPileUp.reads),
                    }
                else
                    .
                end
            ) |
            {
                contigs: $trueContigs,
                correctReads: (.computedReads - (.computedReads - $trueReads)),
                computedReads: .computedReads,
                trueReads: $trueReads,
            }
        ) |
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
        jq --raw-output 'select(has("numReferenceContigs")) | .numReferenceContigs')"

    for (( I = 1; I <= NUM_RESULT_CONTIGS; I++ ));
    do
        result_contig_properly_aligns_to_reference "$I" || return 1
    done
}


#-----------------------------------------------------------------------------
# run main (keep at end of file)
#-----------------------------------------------------------------------------

main
