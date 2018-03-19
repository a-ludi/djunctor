/**
    Defines bindinds and utilities to/for the dazzler commands.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.dazzler;

private {
	auto hiddenDbFileSuffixes = [".bps", ".hdr", ".idx"];
}

/**
	Return a list o hidden files associated to every `.dam`/`.db` file. These
	files contain the actual data used in all the computation. Thus, we
	carefully check for their existence.
*/
auto getHiddenDbFiles(string dbFile)
{
	import std.conv;
	import std.algorithm : map;
	import std.path : baseName, chainPath, dirName, withExtension;

    return hiddenDbFileSuffixes
    	.map!(delegate (suffix) {
    		return chainPath(
    			dbFile.dirName,
    			"." ~ dbFile.baseName.withExtension(suffix).to!string,
			);
		});
}
