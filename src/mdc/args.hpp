#ifndef INCLUDE_MDC_ARGS
#define INCLUDE_MDC_ARGS

#include <cmdline.hpp>
#include "mdc/mdc.hpp"

cmdline::parser *get_parser(int argc, char *argv[]) {
  cmdline::parser *args = new cmdline::parser();

  args->set_program_name("MDC");

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

  args->add(
    "fscore",
    'f',
    "outputs f-score instead of confusion matrix"
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
    false
  );

  args->add<string>(
    "path",
    'p',
    "the path in which to look for datasets",
    false
  );

  args->add<string>(
    "training",
    '\0',
    "the name of the traning set",
    false
  );

  args->add<string>(
    "testing",
    '\0',
    "the name of the testing set",
    false
  );

  args->add<double>(
    "learning_rate",
    '\0',
    "the learning rate used to the SGD on the theta exponents",
    false,
    MDC_DEFAULT_ETA
  );

  args->add<double>(
    "momentum",
    '\0',
    "the momentum used to the SGD on the theta exponents",
    false,
    MDC_DEFAULT_ALPHA
  );

  args->add<double>(
    "delta",
    '\0',
    "the exponent of the prototype distance used to separate classes",
    false,
    MDC_DEFAULT_DELTA
  );

  args->add<double>(
    "beta",
    '\0',
    "the impact applied to the class that is the most distant from its prototype",
    false,
    MDC_DEFAULT_BETA
  );

  args->add<double>(
    "omega",
    '\0',
    "the default description length assumed when there is no clue about the data being assessed",
    false,
    MDC_DEFAULT_OMEGA
  );

  args->add<double>(
    "forgetting_factor",
    '\0',
    "the factor by which the samples in the mixture are considered outdated (1 = no forgetting)",
    false,
    MDC_DEFAULT_F
  );

  args->add<double>(
    "sigma",
    '\0',
    "the standard deviation to the noise applied when the covariance matrix is getting close to be singular",
    false,
    MDC_DEFAULT_SIGMA
  );

  args->parse_check(argc, argv);

  return args;
}

#endif
