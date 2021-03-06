#ifndef INCLUDE_GMDL_ARGS
#define INCLUDE_GMDL_ARGS

#include <cmdline.hpp>
#include "gmdl/gmdl.hpp"

cmdline::parser *get_parser(int argc, char *argv[]) {
  cmdline::parser *args = new cmdline::parser();

  args->set_program_name("GMDL");

  args->add(
    "inline",
    'i',
    "pass configurations parameters inline instead of reading json file"
  );

  args->add(
    "stdin",
    '\0',
    "reads from stdin (prioritized)"
  );

  args->add(
    "quiet",
    'q',
    "omits logging when classification fails"
  );

  args->add(
    "online",
    '\0',
    "makes algorithm online"
  );

  args->add(
    "confusion-matrix",
    '\0',
    "prints confusion matrix"
  );

  args->add<int>(
    "label",
    '\0',
    "the index of the column that should be considered the label",
    false
  );

  args->add<string>(
    "labels",
    '\0',
    "the labels of the dataset comma-separated",
    false
  );

  args->add<string>(
    "set",
    '\0',
    "the key to the set in the config file that should be read",
    false,
    ""
  );

  args->add<string>(
    "path",
    'p',
    "the path in which to look for datasets",
    false,
    ""
  );

  args->add<string>(
    "training",
    '\0',
    "the name of the traning set",
    false,
    ""
  );

  args->add<string>(
    "testing",
    '\0',
    "the name of the testing set",
    false,
    ""
  );

  args->add<string>(
    "config",
    '\0',
    "the path to the config json file",
    false,
    "./config.json"
  );

  args->add<double>(
    "learning_rate",
    '\0',
    "the learning rate used to the SGD on the theta exponents",
    false,
    GMDL_DEFAULT_ETA
  );

  args->add<double>(
    "momentum",
    '\0',
    "the momentum used to the SGD on the theta exponents",
    false,
    GMDL_DEFAULT_ALPHA
  );

  args->add<double>(
    "tau",
    '\0',
    "the exponent of the prototype distance used to separate classes",
    false,
    GMDL_DEFAULT_TAU
  );

  args->add<double>(
    "omega",
    '\0',
    "the default description length assumed when there is no clue about the data being assessed",
    false,
    GMDL_DEFAULT_OMEGA
  );

  args->add<double>(
    "forgetting_factor",
    '\0',
    "the factor by which the samples in the mixture are considered outdated (1 = no forgetting)",
    false,
    GMDL_DEFAULT_F
  );

  args->add<double>(
    "sigma",
    '\0',
    "the standard deviation to the noise applied when the covariance matrix is getting close to be singular",
    false,
    GMDL_DEFAULT_SIGMA
  );

  args->add<int>(
    "dimension",
    '\0',
    "the number of attributes in the dataset (mandatory for online learning)",
    false,
    0
  );

  args->parse_check(argc, argv);

  return args;
}

#endif
