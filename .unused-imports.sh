#!/bin/bash

function main()
{
    for FILE in "$@";
    do
        local UNUSED_SYMBOLS=()

        for SYMBOL in $(find_imported_symbols < "$FILE"); do
            local COUNT=$(one_expression_per_line < "$FILE" | \
                          grep -vE '^\s*import\s+[a-zA-Z0-9_.]+\s*:' | \
                          grep -cE "\b$SYMBOL\b")

            if (( COUNT == 0 )); then
                UNUSED_SYMBOLS+=( $SYMBOL )
            fi
        done

        if (( ${#UNUSED_SYMBOLS[*]} > 0 ));
        then
            local OLD_IFS="$IFS"; IFS="|"
            echo "$FILE: \b(${UNUSED_SYMBOLS[*]})\b"
            IFS="$OLD_IFS"
        fi
    done
}

function find_imported_symbols()
{
    one_expression_per_line | \
    grep -E '^\s*import\s+[a-zA-Z0-9_.]+\s*:' | \
    tr ',;' '  ' | \
    sed -r '{
        s/^.*://
        s/=\s*\S+//g
    }'
}

function one_expression_per_line()
{
    tr '\n' ' ' | \
    tr ';' '\n' | \
    sed -r 's/\bimport\b/\nimport/g'
}

main "$@"
