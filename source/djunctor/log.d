/**
    Central logging facility for djunctor.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.log;

import std.array;
import std.datetime;
import std.format;
import std.stdio;
import core.thread;

private {
    LogLevel minLevel = LogLevel.info;
}

/// Sets the minimum log level to be printed.
void setLogLevel(LogLevel level) nothrow
{
    minLevel = level;
}

LogLevel getLogLevel()
{
    return minLevel;
}

bool shouldLog(LogLevel level)
{
    return level >= minLevel;
}

/**
    Logs a message.
    Params:
        level = The log level for the logged message
        fmt = See http://dlang.org/phobos/std_format.html#format-string
*/
void logDebug(T...)(string fmt, lazy T args) nothrow { log(LogLevel.debug_, fmt, args); }
/// ditto
void logDiagnostic(T...)(string fmt, lazy T args) nothrow { log(LogLevel.diagnostic, fmt, args); }
/// ditto
void logInfo(T...)(string fmt, lazy T args) nothrow { log(LogLevel.info, fmt, args); }
/// ditto
void logWarn(T...)(string fmt, lazy T args) nothrow { log(LogLevel.warn, fmt, args); }
/// ditto
void logError(T...)(string fmt, lazy T args) nothrow { log(LogLevel.error, fmt, args); }

/// ditto
void log(T...)(LogLevel level, string fmt, lazy T args)
nothrow {
    import std.range: chain;

    if( level < minLevel ) return;
    string pref;
    final switch( level ){
        case LogLevel.debug_: pref = "TRACE"; break;
        case LogLevel.diagnostic: pref = "DEBUG"; break;
        case LogLevel.info: pref = "INFO"; break;
        case LogLevel.warn: pref = "WARN"; break;
        case LogLevel.error: pref = "ERROR"; break;
        case LogLevel.fatal: pref = "FATAL"; break;
        case LogLevel.none: assert(false);
    }
    auto threadid = () @trusted { return cast(ulong)cast(void*)Thread.getThis(); } ();
    threadid ^= threadid >> 32;

    try {
        auto txt = appender!string();
        txt.reserve(256);
        formattedWrite(txt, fmt, args);

        if (level >= minLevel) {
            File output;
            if (level == LogLevel.info) () @trusted { output = stdout; } ();
            else () @trusted { output = stderr; } ();
            if (output.isOpen) {
                output.writeln(txt.data);
                output.flush();
            }
        }
    } catch( Exception e ){
        // this is bad but what can we do..
        debug assert(false, e.msg);
    }
}

/// Specifies the log level for a particular log message.
enum LogLevel {
    debug_,
    diagnostic,
    info,
    warn,
    error,
    fatal,
    none
}
