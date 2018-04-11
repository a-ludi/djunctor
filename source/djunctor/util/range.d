/**
    Some additional range functions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.util.range;

import std.range : ElementType, isInputRange;

/**
    This range iterates over fixed-sized chunks of size chunkSize of a source
    range. Source must be an input range. chunkSize must be greater than zero.

    See Also: `std.range.chunks`
    Returns: Range of chunks, ie. `ElementType!Source[]`.
*/
auto arrayChunks(Source)(Source range, in size_t chunkSize) if (isInputRange!Source)
{
    alias Element = ElementType!Source;

    static struct ArrayChunks
    {
        private Source range;
        private const size_t chunkSize;
        private Element[] chunk;

        this(Source range, in size_t chunkSize)
        {
            this.range = range;
            this.chunkSize = chunkSize;
            this.popFront();
        }

        void popFront()
        {
            if (range.empty)
            {
                chunk = null;

                return;
            }

            chunk = new Element[chunkSize];

            foreach (i; 0 .. chunkSize)
            {
                chunk[i] = range.front;

                if (range.empty)
                {
                    chunk = chunk[0 .. i + 1];
                    break;
                }
                else
                {
                    range.popFront();
                }
            }
        }

        @property Element[] front()
        {
            return chunk;
        }

        @property bool empty()
        {
            return chunk is null;
        }
    }

    assert(chunkSize > 0, "chunkSize must be greater than zero");

    return ArrayChunks(range, chunkSize);
}

///
unittest
{
    import std.array : array;
    import std.range : iota;

    auto chunks = iota(10).arrayChunks(2);
    assert(chunks.array == [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9]]);
}
