/**
    Application entry point.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module app;

import std.conv;
import std.stdio;


version (Posix)
{
    version (TestGenerator)
    {
        /// Start `testgen` with the given set of arguments.
        int main(string[] args)
        {
            import testgen.commandline;

            return runTestGenCommandline(args);
        }
    }

    version (Djunctor)
    {
        /// Start `djunctor` with the given set of arguments.
        int main(string[] args)
        {
            import djunctor.commandline;

            return runDjunctorCommandline(args);
        }
    }
}
else
{
    static assert(0, "not compatible with non-POSIX systems.");
}
