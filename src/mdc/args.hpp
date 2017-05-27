#ifndef INCLUDE_MDC_ARGS
#define INCLUDE_MDC_ARGS

#include <cmdline.hpp>

cmdline::parser *get_parser(int argc, char *argv[]) {
  cmdline::parser *args = new cmdline::parser();

  args->add(
    "inline",
    'i',
    "pass configurations parameters inline instead of reading json file"
  );

  args->add(
    "quiet",
    'q',
    "omits logging when classification fails"
  );

  args->parse_check(argc, argv);

  return args;
}

#endif
