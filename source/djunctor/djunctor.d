/**
    This is the main algorithm of this package.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.djunctor;

import djunctor.commandline : Options;
import djunctor.log;
import std.conv;
import std.typecons : Tuple;

template AlignmentContainer(R)
{
    alias AlignmentContainer = Tuple!(R, "a2b", R, "b2a",);
}

/// Start the `djunctor` alorithm with preprocessed options.
void runWithOptions(in ref Options options)
{
    import djunctor.dazzler : provideDamFileInWorkdir, setWorkdir;

    logInfo("starting");

    setWorkdir(options.workdir);
    provideDamFileInWorkdir(options.refFile, options.provideMethod);
    provideDamFileInWorkdir(options.readsFile, options.provideMethod);
}
